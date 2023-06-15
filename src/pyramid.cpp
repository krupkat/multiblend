// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mb/pyramid.h"

#include <cmath>
#include <cstring>

#include <simde/x86/sse4.1.h>

#include "mb/aligned_ptr.h"
#include "mb/pnger.h"
#include "mb/threadpool.h"

namespace multiblend {

/***********************************************************************
 * Constructor/destructor
 ***********************************************************************/
Pyramid::Pyramid(int width, int height, int _levels, int x, int y,
                 mt::ThreadpoolPtr threadpool)
    : threadpool_{threadpool} {
  int n_levels = DefaultNumLevels(width, height);
  if (_levels != 0) {
    if (_levels < 0) {
      n_levels -= _levels;
    } else {
      n_levels = _levels;
    }
  }

  int b;
  int req_alignment;

  b = 0;
  req_alignment = 2;

  for (int n = 0; n < n_levels; ++n) {
    bool x_shift = ((x - b) & (req_alignment - 1)) != 0;
    bool y_shift = ((y - b) & (req_alignment - 1)) != 0;
    int pitch = (width + static_cast<int>(x_shift) + 7) & ~7;
    std::size_t bytes = (std::size_t)pitch *
                        ((height + static_cast<int>(y_shift) + 3) & ~3) *
                        sizeof(float);
    total_bytes_ += bytes;

    levels_.push_back(Level{n,
                            width,
                            height,
                            pitch,
                            pitch >> 2,
                            bytes,
                            {nullptr},
                            x,
                            y,
                            x_shift,
                            y_shift,
                            n != 0 ? levels_[n - 1].x_shift : false,
                            n != 0 ? levels_[n - 1].m128_pitch : 0});

    x -= (static_cast<int>(x_shift) << n);
    y -= (static_cast<int>(y_shift) << n);

    x -= req_alignment;
    y -= req_alignment;

    b -= req_alignment;

    req_alignment <<= 1;

    width = (width + static_cast<int>(x_shift) + 6) >> 1;
    height = (height + static_cast<int>(y_shift) + 6) >> 1;
  }

  for (auto level = levels_.begin(); level < levels_.end(); ++level) {
    level->bands.push_back(0);
    if (level->height >
        threadpool_->GetNThreads() * 4) {  // bands should be at least four
                                           // pixels high to be safe for shrink
      for (int i = 1; i < threadpool_->GetNThreads(); ++i) {
        int b = ((int)(level->height * i * (1.0 / threadpool_->GetNThreads())) &
                 ~3);  // was &~1 to make all bands even height (not sure of
                       // reason); changed to &~3 for vertical SSE processing
                       // (see fastblur)
        level->bands.push_back(b);
      }
    }

    level->bands.push_back(level->height);
  }

  lines_.resize(threadpool_->GetNThreads());

  for (int i = 0; i < threadpool_->GetNThreads(); ++i) {
    lines_[i] = memory::AllocAlignedM128(levels_[0].pitch * sizeof(float));
  }
}

/***********************************************************************
 * copiers
 ***********************************************************************/
void Pyramid::set_lut(int bits, bool gamma) {
  if (lut_bits_ < bits || lut_gamma_ != gamma || (lut_.empty())) {
    lut_.resize(1 << bits);
    lut_bits_ = bits;
    lut_gamma_ = gamma;
  }

  unsigned int l = 1 << bits;
  if (gamma) {
    for (unsigned int i = 0; i < l; ++i) {
      lut_[i] = (float)(i * i);
    }
  } else {
    for (unsigned int i = 0; i < l; ++i) {
      lut_[i] = (float)i;
    }
  }
}

void Pyramid::Copy(uint8_t* src_p, int step, int pitch, bool gamma, int bits) {
  if (step > 1) {
    set_lut(bits, gamma);

    auto tasks = mt::MultiFuture<void>{};
    for (int t = 0; t < (int)levels_[0].bands.size() - 1; ++t) {
      switch (bits) {
        case 8:
          tasks.push_back(threadpool_->Queue([=, this] {
            CopyInterleavedThread_8bit(src_p, step, pitch, levels_[0].bands[t],
                                       levels_[0].bands[t + 1]);
          }));
          break;
        case 16:
          tasks.push_back(threadpool_->Queue([=, this] {
            CopyInterleavedThread_16bit((uint16_t*)src_p, step, pitch,
                                        levels_[0].bands[t],
                                        levels_[0].bands[t + 1]);
          }));
          break;
        case 32:
          break;
      }
    }
    tasks.wait();
    tasks.get();
  } else {
    auto tasks = mt::MultiFuture<void>{};
    for (int t = 0; t < (int)levels_[0].bands.size() - 1; ++t) {
      // planar (only slight improvement with MT, but increases CPU usage)
      switch (bits) {
        case 8:
          tasks.push_back(threadpool_->Queue([=, this] {
            CopyPlanarThread_8bit(src_p, pitch, gamma, levels_[0].bands[t],
                                  levels_[0].bands[t + 1]);
          }));
          break;
        case 10:
        case 12:
        case 14:
        case 16:
          tasks.push_back(threadpool_->Queue([=, this] {
            CopyPlanarThread_16bit((uint16_t*)src_p, pitch, gamma,
                                   levels_[0].bands[t],
                                   levels_[0].bands[t + 1]);
          }));
          break;
        case 32:
          tasks.push_back(threadpool_->Queue([=, this] {
            CopyPlanarThread_32bit((simde__m128*)src_p, pitch, gamma,
                                   levels_[0].bands[t],
                                   levels_[0].bands[t + 1]);
          }));
          break;
      }
    }
    tasks.wait();
    tasks.get();
  }
}

void Pyramid::CopyInterleavedThread_8bit(uint8_t* src_p, int step, int pitch,
                                         int sy, int ey) {
  int x;
  int y;
  uint8_t* src_pp;

  src_p += static_cast<ptrdiff_t>(sy) * pitch;

  float* p_p = levels_[0].data.get();
  p_p += static_cast<ptrdiff_t>(sy) * levels_[0].pitch;

  for (y = sy; y < ey; ++y) {
    src_pp = src_p;
    for (x = 0; x < levels_[0].width; ++x) {
      p_p[x] = lut_[*src_pp];
      src_pp += step;
    }
    for (; x < levels_[0].pitch; ++x) {
      p_p[x] = p_p[x - 1];  // this was commented out, not sure if required
    }
    p_p += levels_[0].pitch;
    src_p += pitch;
  }
}

void Pyramid::CopyInterleavedThread_16bit(uint16_t* src_p, int step, int pitch,
                                          int sy, int ey) {
  int x;
  int y;
  uint16_t* src_pp;

  src_p += static_cast<ptrdiff_t>(sy) * pitch;

  float* p_p = levels_[0].data.get();
  p_p += static_cast<ptrdiff_t>(sy) * levels_[0].pitch;

  for (y = sy; y < ey; ++y) {
    src_pp = src_p;
    for (x = 0; x < levels_[0].width; ++x) {
      p_p[x] = lut_[*src_pp];
      src_pp += step;
    }
    for (; x < levels_[0].pitch; ++x) {
      p_p[x] = p_p[x - 1];
    }
    p_p += levels_[0].pitch;
    src_p += pitch;
  }
}

void Pyramid::CopyPlanarThread_8bit(uint8_t* src_p, int pitch, bool gamma,
                                    int sy, int ey) {
  int x;
  int y;
  simde__m128i pixels;
  simde__m128 fpixels;
  simde__m128i shuffle =
      simde_mm_set_epi32(0x80808003, 0x80808002, 0x80808001, 0x80808000);
  auto* rp_p = (simde__m128*)levels_[0].data.get();
  simde__m128* p_p;
  simde__m128i* src_pp_m;
  int* src_pp_i;
  uint8_t* src_pp_b;

  int sixteens = levels_[0].width >> 4;
  int fours = (levels_[0].width & 0xf) >> 2;
  int ones = (levels_[0].width & 3);
  int extras = levels_[0].pitch - levels_[0].width;

  src_p += static_cast<ptrdiff_t>(sy) * pitch;
  rp_p += static_cast<ptrdiff_t>(sy) * levels_[0].m128_pitch;

  int g;

  for (y = sy; y < ey; ++y) {
    src_pp_m = (simde__m128i*)src_p;
    p_p = rp_p;
    if (gamma) {
      for (x = 0; x < sixteens; ++x) {
        pixels = simde_mm_loadu_si128(src_pp_m++);
        fpixels = simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle));
        simde_mm_store_ps((float*)p_p++, simde_mm_mul_ps(fpixels, fpixels));

        pixels = simde_mm_srli_si128(pixels, 4);
        fpixels = simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle));
        simde_mm_store_ps((float*)p_p++, simde_mm_mul_ps(fpixels, fpixels));

        pixels = simde_mm_srli_si128(pixels, 4);
        fpixels = simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle));
        simde_mm_store_ps((float*)p_p++, simde_mm_mul_ps(fpixels, fpixels));

        pixels = simde_mm_srli_si128(pixels, 4);
        fpixels = simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle));
        simde_mm_store_ps((float*)p_p++, simde_mm_mul_ps(fpixels, fpixels));
      }
    } else {
      for (x = 0; x < sixteens; ++x) {
        pixels = simde_mm_loadu_si128(src_pp_m++);
        simde_mm_store_ps(
            (float*)p_p++,
            simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle)));

        pixels = simde_mm_srli_si128(pixels, 4);
        simde_mm_store_ps(
            (float*)p_p++,
            simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle)));

        pixels = simde_mm_srli_si128(pixels, 4);
        simde_mm_store_ps(
            (float*)p_p++,
            simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle)));

        pixels = simde_mm_srli_si128(pixels, 4);
        simde_mm_store_ps(
            (float*)p_p++,
            simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle)));
      }
    }

    src_pp_i = (int*)src_pp_m;

    for (x = 0; x < fours; ++x) {
      if (gamma) {
        fpixels = simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(
            simde_mm_cvtsi32_si128(*src_pp_i++), shuffle));
        simde_mm_store_ps((float*)p_p++, simde_mm_mul_ps(fpixels, fpixels));
      } else {
        simde_mm_store_ps((float*)p_p++,
                          simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(
                              simde_mm_cvtsi32_si128(*src_pp_i++), shuffle)));
      }
    }

    src_pp_b = (uint8_t*)src_pp_i;
    auto* fp_p = (float*)p_p;

    for (x = 0; x < ones; ++x) {
      if (gamma) {
        g = *src_pp_b++;
        *fp_p++ = (float)g * g;
      } else {
        *fp_p++ = (float)*src_pp_b++;
      }
    }

    float extra = fp_p[-1];
    for (x = 0; x < extras; ++x) {
      *fp_p++ = extra;
    }

    src_p += pitch;
    rp_p += levels_[0].m128_pitch;
  }
}

void Pyramid::CopyPlanarThread_16bit(uint16_t* src_p, int pitch, bool gamma,
                                     int sy, int ey) {
  int x;
  int y;

  simde__m128i pixels;
  simde__m128 fpixels;
  simde__m128i shuffle1 =
      simde_mm_set_epi32(0x80800706, 0x80800504, 0x80800302, 0x80800100);
  simde__m128i shuffle2 =
      simde_mm_set_epi32(0x80800f0e, 0x80800d0c, 0x80800b0a, 0x80800908);
  auto* rp_p = (simde__m128*)levels_[0].data.get();
  simde__m128* p_p;
  simde__m128i* src_pp_m;
  uint16_t* src_pp_w;

  int eights = levels_[0].width >> 3;
  int ones = (levels_[0].width & 7);
  int extras = levels_[0].pitch - levels_[0].width;

  src_p += static_cast<ptrdiff_t>(sy) * pitch;
  rp_p += static_cast<ptrdiff_t>(sy) * levels_[0].m128_pitch;

  int g;

  for (y = sy; y < ey; ++y) {
    src_pp_m = (simde__m128i*)src_p;
    p_p = rp_p;
    if (gamma) {
      for (x = 0; x < eights; ++x) {
        pixels = simde_mm_loadu_si128(src_pp_m++);
        fpixels = simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle1));
        simde_mm_store_ps((float*)p_p++, simde_mm_mul_ps(fpixels, fpixels));

        fpixels = simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle2));
        simde_mm_store_ps((float*)p_p++, simde_mm_mul_ps(fpixels, fpixels));
      }
    } else {
      for (x = 0; x < eights; ++x) {
        pixels = simde_mm_loadu_si128(src_pp_m++);
        simde_mm_store_ps(
            (float*)p_p++,
            simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle1)));

        simde_mm_store_ps(
            (float*)p_p++,
            simde_mm_cvtepi32_ps(simde_mm_shuffle_epi8(pixels, shuffle2)));
      }
    }

    src_pp_w = (uint16_t*)src_pp_m;
    auto* fp_p = (float*)p_p;

    for (x = 0; x < ones; ++x) {
      if (gamma) {
        g = *src_pp_w++;
        *fp_p++ = (float)g * g;
      } else {
        *fp_p++ = (float)*src_pp_w++;
      }
    }

    float extra = fp_p[-1];
    for (x = 0; x < extras; ++x) {
      *fp_p++ = extra;
    }

    src_p += pitch;
    rp_p += levels_[0].m128_pitch;
  }
}

void Pyramid::CopyPlanarThread_32bit(simde__m128* src_p, int pitch, bool gamma,
                                     int sy, int ey) {
  int x;
  int y;

  pitch >>= 2;

  auto* rp_p = (simde__m128*)levels_[0].data.get();

  src_p += static_cast<ptrdiff_t>(sy) * pitch;
  rp_p += static_cast<ptrdiff_t>(sy) * levels_[0].m128_pitch;

  if (gamma) {
    int fours = (levels_[0].width + 3) >> 2;

    for (y = sy; y < ey; ++y) {
      for (x = 0; x < fours; ++x) {
        simde__m128 pixels = simde_mm_load_ps((float*)&src_p[x]);
        simde_mm_store_ps((float*)&rp_p[x], simde_mm_mul_ps(pixels, pixels));
      }
      float copy = ((float*)rp_p)[levels_[0].width - 1];
      for (x = levels_[0].width; x < levels_[0].pitch; ++x) {
        ((float*)rp_p)[x] = copy;
      }
      src_p += pitch;
      rp_p += levels_[0].m128_pitch;
    }
  } else {
    int copy_bytes = levels_[0].width << 2;

    for (y = sy; y < ey; ++y) {
      memcpy(rp_p, src_p, copy_bytes);
      float copy = ((float*)rp_p)[levels_[0].width - 1];
      for (x = levels_[0].width; x < levels_[0].pitch; ++x) {
        ((float*)rp_p)[x] = copy;
      }
      src_p += pitch;
      rp_p += levels_[0].m128_pitch;
    }
  }
}

/***********************************************************************
 * subsample
 ***********************************************************************/
void Pyramid::Subsample(int sub_w, int sub_h, Pyramid* source) {
  int x;
  int y;
  int p = 0;
  auto* in = (simde__m128*)source->levels_[0].data.get();
  auto* out = (simde__m128*)levels_[0].data.get();
  auto& temp_lines = source->lines_;
  simde__m128* line = temp_lines[0].get();
  int m128_pitch_in = source->levels_[0].m128_pitch;
  int m128_pitch_out = levels_[0].m128_pitch;
  int mid_pitch = sub_w == 2 ? ((m128_pitch_in >> 1) + 1) & ~1 : m128_pitch_in;
  simde__m128 three = simde_mm_set_ps1(3);
  simde__m128 four = simde_mm_set_ps1(4);
  simde__m128 mul = simde_mm_set_ps1((sub_h != 0 ? 1.0f / 8 : 1.0f) *
                                     (sub_w == 2 ? 1.0f / 64 : 1.0f / 8));

  for (y = 0; y < levels_[0].height; ++y) {
    if (sub_h != 0) {
      if (y == 0) {
        for (x = 0; x < m128_pitch_in; ++x) {
          simde_mm_store_ps(
              (float*)&temp_lines[0][x],
              simde_mm_add_ps(
                  simde_mm_mul_ps(simde_mm_load_ps((float*)&in[x]), four),
                  simde_mm_add_ps(
                      simde_mm_mul_ps(
                          simde_mm_load_ps((float*)&in[x + m128_pitch_in]),
                          three),
                      simde_mm_load_ps(
                          (float*)&in[x + (m128_pitch_in << 1)]))));
        }
      } else if (y == levels_[0].height - 1) {
        for (x = 0; x < m128_pitch_in; ++x) {
          simde_mm_store_ps(
              (float*)&temp_lines[0][x],
              simde_mm_add_ps(
                  simde_mm_mul_ps(
                      simde_mm_load_ps((float*)&in[x + m128_pitch_in]), four),
                  simde_mm_add_ps(
                      simde_mm_mul_ps(simde_mm_load_ps((float*)&in[x]), three),
                      simde_mm_load_ps((float*)&in[x - m128_pitch_in]))));
        }
      } else {
        for (x = 0; x < m128_pitch_in; ++x) {
          simde_mm_store_ps(
              (float*)&temp_lines[0][x],
              simde_mm_add_ps(
                  simde_mm_add_ps(
                      simde_mm_load_ps((float*)&in[x - m128_pitch_in]),
                      simde_mm_load_ps((float*)&in[x + (m128_pitch_in << 1)])),
                  simde_mm_mul_ps(
                      simde_mm_add_ps(
                          simde_mm_load_ps((float*)&in[x]),
                          simde_mm_load_ps((float*)&in[x + m128_pitch_in])),
                      three)));
        }
      }
      in += m128_pitch_in << 1;
    } else {
      line = (simde__m128*)in;
      in += m128_pitch_in;
    }
    switch (sub_w) {
      case 2:
        Subsample_Squeeze(line, line, m128_pitch_in, mid_pitch, nullptr);
      case 1:
        Subsample_Squeeze(line, out, mid_pitch, m128_pitch_out, &mul);
    }
    out += m128_pitch_out;
  }
}

void Pyramid::Subsample_Squeeze(simde__m128* in, simde__m128* Out,
                                int m128_pitch_in, int m128_pitch_out,
                                simde__m128* mul) {
  int read = 0;

  int x;
  simde__m128 a;
  simde__m128 b;
  simde__m128 c;
  simde__m128 d;
  simde__m128 e;
  simde__m128 f;
  simde__m128 three = simde_mm_set_ps1(3);

  b = simde_mm_load_ps((float*)&in[read++]);
  a = simde_mm_shuffle_ps(b, b, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
  c = simde_mm_load_ps((float*)&in[read++]);
  d = simde_mm_load_ps((float*)&in[read++]);

  for (x = 0; x < m128_pitch_out; ++x) {
    e = simde_mm_shuffle_ps(a, c, SIMDE_MM_SHUFFLE(0, 0, 3, 3));
    f = simde_mm_shuffle_ps(b, d, SIMDE_MM_SHUFFLE(0, 0, 3, 3));
    e = simde_mm_blend_ps(b, e, 9);
    f = simde_mm_blend_ps(c, f, 9);
    e = simde_mm_shuffle_ps(e, e, SIMDE_MM_SHUFFLE(3, 1, 2, 0));
    f = simde_mm_shuffle_ps(f, f, SIMDE_MM_SHUFFLE(3, 1, 2, 0));
    e = simde_mm_hadd_ps(e, f);

    f = simde_mm_hadd_ps(b, c);
    if (mul != nullptr) {
      simde_mm_store_ps(
          (float*)&Out[x],
          simde_mm_mul_ps(simde_mm_add_ps(e, simde_mm_mul_ps(f, three)), *mul));
    } else {
      simde_mm_store_ps((float*)&Out[x],
                        simde_mm_add_ps(e, simde_mm_mul_ps(f, three)));
    }
    a = c;
    b = d;
    if (read < m128_pitch_in - 1) {
      c = simde_mm_load_ps((float*)&in[read++]);
      d = simde_mm_load_ps((float*)&in[read++]);
    } else {
      if (read < m128_pitch_in) {
        c = simde_mm_load_ps((float*)&in[read++]);
        d = simde_mm_shuffle_ps(c, c, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
      } else {
        c = simde_mm_shuffle_ps(b, b, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
        d = c;
      }
    }
  }
}

/***********************************************************************
 * shrink (gaussian)
 ***********************************************************************/
void Pyramid::Shrink() {
  int l;
  simde__m128* hi;
  simde__m128* lo;
  const simde__m128 four = simde_mm_set_ps1(4);
  const simde__m128 six = simde_mm_set_ps1(6);
  const simde__m128 eleven = simde_mm_set_ps1(11);
  const simde__m128 fifteen = simde_mm_set_ps1(15);
  const simde__m128 _16th = simde_mm_set_ps1(1.0 / 16);
  const simde__m128 _256th = simde_mm_set_ps1(1.0 / 256);

  for (l = 0; l < (int)levels_.size() - 1; ++l) {
    hi = (simde__m128*)levels_[l].data.get();
    lo = (simde__m128*)levels_[l + 1].data.get();

    int height_odd =
        (levels_[l].height & 1) ^ static_cast<int>(levels_[l].y_shift);
    int first_bad_line = levels_[l + 1].height - (3 - height_odd);

    auto tasks = mt::MultiFuture<void>{};
    for (int t = 0; t < (int)levels_[l + 1].bands.size() - 1; ++t) {
      tasks.push_back(threadpool_->Queue([=, this] {
        ShrinkThread(lines_[t].get(), hi, lo, levels_[l].m128_pitch,
                     levels_[l + 1].m128_pitch, first_bad_line, height_odd,
                     levels_[l + 1].bands[t], levels_[l + 1].bands[t + 1],
                     levels_[l].x_shift, levels_[l].y_shift);
      }));
    }
    tasks.wait();
    tasks.get();
  }
}

void Pyramid::ShrinkThread(simde__m128* line, simde__m128* hi, simde__m128* lo,
                           int m128_pitch_hi, int m128_pitch_lo,
                           int first_bad_line, int height_odd, int sy, int ey,
                           const bool x_shift, const bool y_shift) {
  int x;
  int y;

  const simde__m128 four = simde_mm_set_ps1(4);
  const simde__m128 six = simde_mm_set_ps1(6);
  const simde__m128 eleven = simde_mm_set_ps1(11);
  const simde__m128 fifteen = simde_mm_set_ps1(15);
  const simde__m128 _16th = simde_mm_set_ps1(1.0 / 16);
  const simde__m128 _256th = simde_mm_set_ps1(1.0 / 256);

  // line 0
  if (sy == 0) {
    memcpy(line, hi, m128_pitch_hi << 4);
    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _16th, x_shift);
    lo += m128_pitch_lo;

    if (!y_shift) {
      // line 1
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] = simde_mm_add_ps(
            simde_mm_add_ps(simde_mm_mul_ps(hi[x], eleven),
                            simde_mm_mul_ps(hi[x + m128_pitch_hi], four)),
            hi[x + (m128_pitch_hi << 1)]);
      }
      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      lo += m128_pitch_lo;
      hi += (m128_pitch_hi << 1);
      sy = 2;
    } else {
      // line 1
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] = simde_mm_add_ps(simde_mm_mul_ps(hi[x], fifteen),
                                  hi[x + m128_pitch_hi]);
      }
      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      lo += m128_pitch_lo;
      hi += m128_pitch_hi;
      // line 2
      for (x = 0; x < m128_pitch_hi; ++x) {
        simde_mm_store_ps(
            (float*)&line[x],
            simde_mm_add_ps(
                simde_mm_add_ps(
                    simde_mm_add_ps(hi[x - m128_pitch_hi],
                                    hi[x + (m128_pitch_hi << 1)]),
                    simde_mm_mul_ps(simde_mm_add_ps(hi[x - m128_pitch_hi],
                                                    hi[x + m128_pitch_hi]),
                                    four)),
                simde_mm_mul_ps(hi[x], six)));
      }
      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      lo += m128_pitch_lo;
      hi += (m128_pitch_hi << 1);
      sy = 3;
    }
  } else {
    lo += static_cast<ptrdiff_t>(m128_pitch_lo) * sy;
    hi +=
        static_cast<ptrdiff_t>(m128_pitch_hi) *
        (((sy - 1) << 1) -
         static_cast<int>(
             y_shift));  // because hi has a missing line compared to lo, lo's
                         // line 0 corresponds to hi's line -2, 1=>0, 2=>2, 3=>4
  }

  ey = (std::min)(first_bad_line, ey);

  // good lines
  for (y = sy; y < ey; ++y) {  // was y < first_bad_line
    for (x = 0; x < m128_pitch_hi; ++x) {
      simde_mm_store_ps(
          (float*)&line[x],
          simde_mm_add_ps(
              simde_mm_add_ps(
                  simde_mm_add_ps(hi[x - (m128_pitch_hi << 1)],
                                  hi[x + (m128_pitch_hi << 1)]),
                  simde_mm_mul_ps(simde_mm_add_ps(hi[x - m128_pitch_hi],
                                                  hi[x + m128_pitch_hi]),
                                  four)),
              simde_mm_mul_ps(hi[x], six)));
    }
    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
    lo += m128_pitch_lo;
    hi += (m128_pitch_hi << 1);
  }

  if (y == first_bad_line) {  // final block
    // prepenultimate line
    if (height_odd == 0) {
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] = simde_mm_add_ps(
            simde_mm_add_ps(
                simde_mm_add_ps(hi[x - (m128_pitch_hi << 1)],
                                hi[x + m128_pitch_hi]),
                simde_mm_mul_ps(simde_mm_add_ps(hi[x - m128_pitch_hi],
                                                hi[x + m128_pitch_hi]),
                                four)),
            simde_mm_mul_ps(hi[x], six));
      }

      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      ++y;
      lo += m128_pitch_lo;
      hi += (m128_pitch_hi << 1);

      // this case moved from block below
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] =
            simde_mm_add_ps(hi[x - (m128_pitch_hi << 1)],
                            simde_mm_mul_ps(hi[x - m128_pitch_hi], fifteen));
      }
    } else {
      // penultimate line
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] = simde_mm_add_ps(
            simde_mm_add_ps(simde_mm_mul_ps(hi[x], eleven),
                            simde_mm_mul_ps(hi[x - m128_pitch_hi], four)),
            hi[x - (m128_pitch_hi << 1)]);
      }
      // other case removed from here, moved into block above
    }

    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
    ++y;
    lo += m128_pitch_lo;
    hi += (m128_pitch_hi << 1);

    // last line
    hi -= static_cast<ptrdiff_t>(m128_pitch_hi) * (3 - height_odd);
    memcpy(line, hi, m128_pitch_hi << 4);
    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _16th, x_shift);
    ++y;
  }
}

void Pyramid::Squeeze(simde__m128* line, simde__m128* lo, int m128_pitch_lo,
                      int m128_pitch_hi, simde__m128 final_mul, bool x_shift) {
  int hi_x = 0;
  int lo_x = 0;

  simde__m128 a;
  simde__m128 b;
  simde__m128 c;
  simde__m128 d;
  simde__m128 e;
  simde__m128 f;
  simde__m128 g;
  simde__m128 h;
  simde__m128 i;
  simde__m128 j;

  const simde__m128 four = simde_mm_set_ps1(4);
  const simde__m128 six = simde_mm_set_ps1(6);

  if (x_shift) {
    memmove(&((float*)line)[1], line, (m128_pitch_hi << 4) - 4);
  }

  while (lo_x < m128_pitch_lo) {
    if (hi_x >= m128_pitch_hi) {  // was >= ... + 1
      b = simde_mm_shuffle_ps(a, a, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
      c = b;
    } else {
      b = line[hi_x++];
      c = line[hi_x++];
      if (lo_x == 0) {
        a = simde_mm_shuffle_ps(b, b, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
      }
    }

    // a = EFGH
    // b = IJKL
    // c = MNOP

    // shuffle to four pairs of outer pixels
    f = simde_mm_shuffle_ps(a, b, SIMDE_MM_SHUFFLE(2, 0, 2, 0));  // EGIK
    g = simde_mm_shuffle_ps(b, c, SIMDE_MM_SHUFFLE(2, 0, 2, 0));  // IKMO
    d = simde_mm_shuffle_ps(f, f, SIMDE_MM_SHUFFLE(3, 1, 2, 0));  // EIGK
    e = simde_mm_shuffle_ps(g, g, SIMDE_MM_SHUFFLE(3, 1, 2, 0));  // IMKO
    d = simde_mm_hadd_ps(d, e);

    // shuffle to four pairs of inner pixels
    h = simde_mm_shuffle_ps(
        a, b,
        SIMDE_MM_SHUFFLE(1, 1, 3, 1));  // FHJJ // prev BDFF (should be BDDF)
    i = simde_mm_shuffle_ps(
        b, c,
        SIMDE_MM_SHUFFLE(1, 1, 3, 1));  // JLNN // prev FHJJ (should be FHHJ)
    h = simde_mm_shuffle_ps(h, h, SIMDE_MM_SHUFFLE(3, 1, 1, 0));  // FHHJ
    i = simde_mm_shuffle_ps(i, i, SIMDE_MM_SHUFFLE(3, 1, 1, 0));  // JLLN
    h = simde_mm_mul_ps(four, simde_mm_hadd_ps(h, i));

    // shuffle to four central pixels
    j = simde_mm_mul_ps(
        six, simde_mm_shuffle_ps(f, g, SIMDE_MM_SHUFFLE(2, 1, 2, 1)));  // GIKM

    // store
    simde_mm_store_ps(
        (float*)&lo[lo_x++],
        simde_mm_mul_ps(simde_mm_add_ps(d, simde_mm_add_ps(h, j)), final_mul));
    a = c;
  }
}

/***********************************************************************
 * laplace (gaussian)
 ***********************************************************************/
void Pyramid::LaplaceCollapse(int n_levels, bool Collapse) {
  int j;
  int l;

  for (j = 0; j < n_levels - 1; ++j) {
    if (Collapse) {
      l = (n_levels - 2) - j;
    } else {
      l = j;
    }

    auto tasks = mt::MultiFuture<void>{};
    for (int t = 0; t < (int)levels_[l].bands.size() - 1; ++t) {
      tasks.push_back(threadpool_->Queue([=, this] {
        LaplaceThreadWrapper(&levels_[l], &levels_[l + 1], levels_[l].bands[t],
                             levels_[l].bands[t + 1]);
      }));
    }
    tasks.wait();
    tasks.get();
  }
}

void Pyramid::LaplaceExpand(simde__m128* hi, simde__m128* lo, int m128_pitch_hi,
                            int m128_pitch_lo) {
  simde__m128 p;
  simde__m128 q;

  const simde__m128 expand0 = simde_mm_set_ps(0, 0.125f, 0.75f, 0.125f);
  const simde__m128 expand1 = simde_mm_set_ps(0, 0.5f, 0.5f, 0);
  const simde__m128 expand2 = simde_mm_set_ps(0.125, 0.75f, 0.125f, 0);
  const simde__m128 expand3 = simde_mm_set_ps(0.5f, 0.5f, 0, 0);

  int x_hi = 0;
  int x_lo = 0;

  p = simde_mm_load_ps((float*)lo);
  ++x_lo;

  while (x_hi < m128_pitch_hi) {
    simde_mm_store_ps(
        (float*)&hi[x_hi],
        simde_mm_hadd_ps(simde_mm_hadd_ps(simde_mm_mul_ps(p, expand0),
                                          simde_mm_mul_ps(p, expand1)),
                         simde_mm_hadd_ps(simde_mm_mul_ps(p, expand2),
                                          simde_mm_mul_ps(p, expand3))));

    ++x_hi;

    if (x_lo < m128_pitch_lo) {
      q = lo[x_lo];  // simde_mm_load_ps((float*)&lo[x_lo]);
    } else if (x_lo == m128_pitch_lo) {
      q = simde_mm_shuffle_ps(q, q, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
    }
    ++x_lo;
    p = simde_mm_shuffle_ps(p, q, SIMDE_MM_SHUFFLE(1, 0, 3, 2));

    simde_mm_store_ps(
        (float*)&hi[x_hi],
        simde_mm_hadd_ps(simde_mm_hadd_ps(simde_mm_mul_ps(p, expand0),
                                          simde_mm_mul_ps(p, expand1)),
                         simde_mm_hadd_ps(simde_mm_mul_ps(p, expand2),
                                          simde_mm_mul_ps(p, expand3))));

    p = q;

    ++x_hi;
  }
}

void Pyramid::LaplaceExpandShifted(simde__m128* hi, simde__m128* lo,
                                   int m128_pitch_hi, int m128_pitch_lo) {
  simde__m128 p;
  simde__m128 q;
  simde__m128 t;

  const simde__m128 expand0 = simde_mm_set_ps(0, 0.5f, 0.5f, 0);
  const simde__m128 expand1 = simde_mm_set_ps(0.125, 0.75f, 0.125f, 0);
  const simde__m128 expand2 = simde_mm_set_ps(0.5f, 0.5f, 0, 0);
  const simde__m128 expand3 = simde_mm_set_ps(0, 0.125f, 0.75f, 0.125f);

  int x_hi = 0;
  int x_lo = 0;

  p = simde_mm_load_ps((float*)lo);
  ++x_lo;

  while (x_hi < m128_pitch_hi) {
    t = p;

    if (x_lo < m128_pitch_lo) {
      q = lo[x_lo];
    } else if (x_lo == m128_pitch_lo) {
      q = simde_mm_shuffle_ps(q, q, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
    }
    ++x_lo;
    p = simde_mm_shuffle_ps(p, q, SIMDE_MM_SHUFFLE(1, 0, 3, 2));

    simde_mm_store_ps(
        (float*)&hi[x_hi],
        simde_mm_hadd_ps(simde_mm_hadd_ps(simde_mm_mul_ps(t, expand0),
                                          simde_mm_mul_ps(t, expand1)),
                         simde_mm_hadd_ps(simde_mm_mul_ps(t, expand2),
                                          simde_mm_mul_ps(p, expand3))));

    ++x_hi;

    simde_mm_store_ps(
        (float*)&hi[x_hi],
        simde_mm_hadd_ps(simde_mm_hadd_ps(simde_mm_mul_ps(p, expand0),
                                          simde_mm_mul_ps(p, expand1)),
                         simde_mm_hadd_ps(simde_mm_mul_ps(p, expand2),
                                          simde_mm_mul_ps(q, expand3))));

    p = q;

    ++x_hi;
  }
}

void Pyramid::LaplaceLine2(simde__m128* hi, simde__m128* temp1,
                           simde__m128* temp2, int m128_pitch) {
  const simde__m128 half = simde_mm_set_ps1(0.5f);

  for (int x = 0; x < m128_pitch; ++x) {
    simde_mm_store_ps(
        (float*)&hi[x],
        simde_mm_sub_ps(
            simde_mm_mul_ps(simde_mm_add_ps(temp1[x], temp2[x]), half), hi[x]));
  }
}

void Pyramid::LaplaceLine3(simde__m128* hi, simde__m128* temp1,
                           simde__m128* temp2, simde__m128* temp3,
                           int m128_pitch) {
  static const simde__m128 eighth = simde_mm_set_ps1(0.125f);
  static const simde__m128 three_quarters = simde_mm_set_ps1(0.75f);

  for (int x = 0; x < m128_pitch; ++x) {
    simde_mm_store_ps(
        (float*)&hi[x],
        simde_mm_sub_ps(
            simde_mm_add_ps(
                simde_mm_mul_ps(simde_mm_add_ps(temp1[x], temp3[x]), eighth),
                simde_mm_mul_ps(temp2[x], three_quarters)),
            hi[x]));
  }
}

simde__m128* GetLine(const Pyramid::Level& level, int y) {
  return (simde__m128*)(level.data.get() +
                        static_cast<ptrdiff_t>(y) * level.pitch);
}

void GetExpandedLine(const Pyramid::Level& level, simde__m128* temp, int y) {
  if (level.upper_x_shift) {
    Pyramid::LaplaceExpandShifted(temp, GetLine(level, y),
                                  level.upper_m128_pitch, level.m128_pitch);
  } else {
    Pyramid::LaplaceExpand(temp, GetLine(level, y), level.upper_m128_pitch,
                           level.m128_pitch);
  }
}

void Pyramid::LaplaceThreadWrapper(Level* upper_level, Level* lower_level,
                                   int sy, int ey) {
  int temp = upper_level->m128_pitch << 4;

  auto temp1 = memory::AllocAlignedM128(temp);
  auto temp2 = memory::AllocAlignedM128(temp);
  auto temp3 = memory::AllocAlignedM128(temp);

  LaplaceThread(upper_level, lower_level, sy, ey, temp1.get(), temp2.get(),
                temp3.get());
}

void Pyramid::LaplaceThread(Level* upper_level, Level* lower_level, int sy,
                            int ey, simde__m128* temp1, simde__m128* temp2,
                            simde__m128* temp3) {
  simde__m128* hi = (simde__m128*)upper_level->data.get() +
                    static_cast<ptrdiff_t>(sy) * upper_level->m128_pitch;

  int lo_y = sy >> 1;

  if (upper_level->y_shift) {
    ++lo_y;
    GetExpandedLine(*lower_level, temp2, lo_y++);
    GetExpandedLine(*lower_level, temp3, lo_y++);
  } else {
    GetExpandedLine(*lower_level, temp1, lo_y++);
    GetExpandedLine(*lower_level, temp2, lo_y++);
  }

  for (int y = sy; y < ey; ++y) {
    if (((y + static_cast<int>(upper_level->y_shift)) & 1) == 0) {
      GetExpandedLine(*lower_level, temp3, lo_y++);

      LaplaceLine3(hi, temp1, temp2, temp3, upper_level->m128_pitch);
    } else {
      // NOLINTNEXTLINE(readability-suspicious-call-argument)
      LaplaceLine2(hi, temp2, temp3, upper_level->m128_pitch);

      simde__m128* temp = temp1;
      temp1 = temp2;
      temp2 = temp3;
      temp3 = temp;
    }

    hi += upper_level->m128_pitch;
  }
}

/***********************************************************************
 * Average (top level)
 ***********************************************************************/
float Pyramid::Average() {
  int x;
  int y;
  int fours = levels_[0].width >> 2;

  simde__m128 m128_total = {0};
  simde__m128 one = simde_mm_set_ps1(1.0f);
  double total = 0;
  double row_total;

  auto* data = (simde__m128*)levels_[0].data.get();

  for (y = 0; y < levels_[0].height; ++y) {
    m128_total = simde_mm_setzero_ps();

    for (x = 0; x < fours; ++x) {
      m128_total = simde_mm_add_ps(m128_total, data[x]);
    }

    m128_total = simde_mm_hadd_ps(m128_total, m128_total);
    m128_total = simde_mm_hadd_ps(m128_total, m128_total);
    row_total = simde_mm_cvtss_f32(m128_total);

    for (x <<= 2; x < levels_[0].width; ++x) {
      row_total += ((float*)data)[x];
    }

    total += row_total;

    data += levels_[0].m128_pitch;
  }

  total /= levels_[0].width;
  total /= levels_[0].height;

  return (float)total;
}

/***********************************************************************
 * Add
 ***********************************************************************/
void Pyramid::Add(float add, int levels) {
  simde__m128 __add = simde_mm_set_ps1(add);

  int lim = (std::min)(levels, (int)levels_.size() - 1);

  for (int l = 0; l < lim; ++l) {
    auto* data = (simde__m128*)levels_[l].data.get();

    auto tasks = mt::MultiFuture<void>{};
    for (int t = 0; t < (int)levels_[l].bands.size() - 1; ++t) {
      tasks.push_back(threadpool_->Queue([=, this]() {
        simde__m128* data =
            (simde__m128*)levels_[l].data.get() +
            static_cast<ptrdiff_t>(levels_[l].bands[t]) * levels_[l].m128_pitch;
        for (int y = levels_[l].bands[t]; y < levels_[l].bands[t + 1]; ++y) {
          for (int x = 0; x < levels_[l].m128_pitch; ++x) {
            simde_mm_store_ps((float*)&data[x],
                              simde_mm_add_ps(__add, data[x]));
          }
          data += levels_[l].m128_pitch;
        }
      }));
    }
    tasks.wait();
    tasks.get();
  }
}

/***********************************************************************
 * Multiply and add
 ***********************************************************************/
void Pyramid::MultiplyAndAdd(float add, float mul, int levels) {
  int i;
  int x;
  int y;
  simde__m128 __add = simde_mm_set_ps1(add);
  simde__m128 __mul = simde_mm_set_ps1(mul);

  int lim = (std::min)(levels, (int)levels_.size() - 1);

  for (i = 0; i < lim; ++i) {
    auto* data = (simde__m128*)levels_[i].data.get();

    for (y = 0; y < levels_[i].height; ++y) {
      for (x = 0; x < levels_[i].m128_pitch; ++x) {
        simde_mm_store_ps(
            (float*)&data[x],
            simde_mm_add_ps(__add, simde_mm_mul_ps(data[x], __mul)));
      }
      data += levels_[i].m128_pitch;
    }
  }
}

/***********************************************************************
 * Multiply, add, clamp
 ***********************************************************************/
void Pyramid::MultiplyAddClamp(float add, float mul, int level) {
  int x;
  int y;
  simde__m128 __add = simde_mm_set_ps1(add);
  simde__m128 __mul = simde_mm_set_ps1(mul);
  simde__m128 __min = simde_mm_set_ps1(1.0f);
  simde__m128 __max = simde_mm_set_ps1(0.0f);

  auto* data = (simde__m128*)levels_[level].data.get();

  for (y = 0; y < levels_[level].height; ++y) {
    for (x = 0; x < levels_[level].m128_pitch; ++x) {
      simde_mm_store_ps(
          (float*)&data[x],
          simde_mm_max_ps(
              simde_mm_min_ps(
                  simde_mm_add_ps(__add, simde_mm_mul_ps(data[x], __mul)),
                  __min),
              __max));
    }
    data += levels_[level].m128_pitch;
  }
}

/***********************************************************************
 * Multiply
 ***********************************************************************/
void Pyramid::Multiply(int level, float mul) {
  if (mul == 1) {
    return;
  }
  if (mul == 0) {
    memset(levels_[level].data.get(), 0,
           static_cast<std::size_t>(levels_[level].height) *
               levels_[level].pitch * sizeof(float));
    return;
  }

  int x;
  int y;
  simde__m128 __mul = simde_mm_set_ps1(mul);

  auto* data = (simde__m128*)levels_[level].data.get();

  for (y = 0; y < levels_[level].height; ++y) {
    for (x = 0; x < levels_[level].m128_pitch; ++x) {
      simde_mm_store_ps(
          (float*)&data[x],
          simde_mm_mul_ps(simde_mm_load_ps((float*)&data[x]), __mul));
    }
    data += levels_[level].m128_pitch;
  }
}

/***********************************************************************
 * multiply_by_pyramid
 ***********************************************************************/
void Pyramid::MultplyByPyramid(Pyramid* b) {
  int x;
  int y;

  for (int l = 0; l < (int)levels_.size() - 1; ++l) {
    auto* data = (simde__m128*)levels_[l].data.get();
    auto* _b = (simde__m128*)b->levels_[l].data.get();

    for (y = 0; y < levels_[l].height; ++y) {
      for (x = 0; x < levels_[l].m128_pitch; ++x) {
        simde_mm_store_ps((float*)&data[x],
                          simde_mm_mul_ps(simde_mm_load_ps((float*)&data[x]),
                                          simde_mm_load_ps((float*)&_b[x])));
      }
      data += levels_[l].m128_pitch;
      _b += levels_[l].m128_pitch;
    }
  }
}

/***********************************************************************
 * blend
 ***********************************************************************/
void Pyramid::Fuse(Pyramid* _b, Pyramid* mask, bool pre = false,
                   int black = 0x00) {
  int l;

  for (l = 0; l < (int)levels_.size(); ++l) {
    //  fuse_thread((simde__m128*)levels[l].data,
    //(simde__m128*)_b->levels[l].data,
    //(simde__m128*)mask->levels[l].data, m128_pitch, 0, levels[l].height, pre,
    // black);

    // fuse doesn't see any gains from multithreading; leave this here as
    // reference

    auto tasks = mt::MultiFuture<void>{};
    for (int t = 0; t < (int)levels_[l].bands.size() - 1; ++t) {
      tasks.push_back(threadpool_->Queue([=, this] {
        FuseThread((simde__m128*)levels_[l].data.get(),
                   (simde__m128*)_b->levels_[l].data.get(),
                   (simde__m128*)mask->levels_[l].data.get(),
                   levels_[l].m128_pitch, levels_[l].bands[t],
                   levels_[l].bands[t + 1], pre, black);
      }));
    }
    tasks.wait();
    tasks.get();
  }
}

void Pyramid::FuseThread(simde__m128* a, simde__m128* b, simde__m128* m,
                         int m128_pitch, int sy, int ey, bool pre, int black) {
  int p;
  int add = sy * m128_pitch;
  int count = (ey - sy) * m128_pitch;

  a += add;
  b += add;
  m += add;

  if (!pre) {
    for (p = 0; p < count; ++p) {
      simde__m128 _a = a[p];
      simde_mm_store_ps(
          (float*)&a[p],
          simde_mm_add_ps(_a,
                          simde_mm_mul_ps(simde_mm_sub_ps(b[p], _a), m[p])));
    }
  } else {
    simde__m128 ones = simde_mm_set_ps1(1.0f);
    simde__m128 blacks = simde_mm_set_ps1((float)black);
    if (black != 0) {
      for (p = 0; p < count; ++p) {
        simde_mm_store_ps(
            (float*)&a[p],
            simde_mm_add_ps(
                blacks,
                simde_mm_add_ps(simde_mm_sub_ps(b[p], blacks),
                                simde_mm_mul_ps(simde_mm_sub_ps(a[p], blacks),
                                                simde_mm_sub_ps(ones, m[p])))));
      }
    } else {
      for (p = 0; p < count; ++p) {
        simde_mm_store_ps(
            (float*)&a[p],
            simde_mm_add_ps(
                b[p], simde_mm_mul_ps(a[p], simde_mm_sub_ps(ones, m[p]))));
      }
    }
  }
}

void Pyramid::Fuse(Pyramid* b, float weight) {
  int l;
  int p;
  simde__m128 w = simde_mm_set_ps1(weight);

  for (l = 0; l < (int)levels_.size(); ++l) {
    auto* _a = (simde__m128*)levels_[l].data.get();
    auto* _b = (simde__m128*)b->levels_[l].data.get();

    int count = levels_[l].height * levels_[l].m128_pitch;
    for (p = 0; p < count; ++p) {
      simde__m128 __a = _a[p];
      simde_mm_store_ps(
          (float*)&_a[p],
          simde_mm_add_ps(__a,
                          simde_mm_mul_ps(simde_mm_sub_ps(_b[p], __a), w)));
    }
  }
}

/***********************************************************************
 * denoise
 ***********************************************************************/
#ifdef PYR_DENOISE
void Pyramid::Denoise(int level, float power, bool gamma) {
  if (power == 0) return;

  int x, y;
  simde__m128 one = simde_mm_set_ps1(1);
  simde__m128 half = simde_mm_set_ps1(0.5f);
  simde__m128 _power = simde_mm_set_ps1(power);
  simde__m128 pi = simde_mm_set_ps1(3.14159265359f);
  simde__m128i andi =
      simde_mm_set_epi32(0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff);
  simde__m128* _and = (simde__m128*)&andi;

  if (gamma) _power = simde_mm_mul_ps(_power, _power);
  _power = simde_mm_div_ps(one, _power);

  simde__m128* data = (simde__m128*)levels[level].data;
  for (y = 0; y < levels[level].height; ++y) {
    for (x = 0; x < levels[level].m128_pitch; ++x) {
      simde__m128 d = data[x];
      simde_mm_store_ps(
          (float*)&data[x],
          simde_mm_mul_ps(
              simde_mm_mul_ps(
                  simde_mm_sub_ps(
                      one,
                      cos_ps(simde_mm_min_ps(
                          simde_mm_and_ps(simde_mm_mul_ps(d, _power), *_and),
                          pi))),
                  half),
              d));
    }
    data += levels[level].m128_pitch;
  }
}
#endif

/***********************************************************************
 * blend (base swap)
 ***********************************************************************/
void Pyramid::Blend(Pyramid* b) {
  if (b->GetNLevels() < GetNLevels()) {
    return;
  }
  memcpy(levels_[GetNLevels() - 1].data.get(),
         b->levels_[GetNLevels() - 1].data.get(),
         static_cast<std::size_t>(levels_[GetNLevels() - 1].height) *
             levels_[GetNLevels() - 1].pitch * sizeof(float));
}

/***********************************************************************
 * approximate gaussian blur
 ***********************************************************************/

void Pyramid::BlurX(float radius, Pyramid* transpose) {
  auto tasks = mt::MultiFuture<void>{};
  for (int i = 0; i < (int)levels_[0].bands.size() - 1; ++i) {
    tasks.push_back(threadpool_->Queue([=, this] {
      BlurXThread(radius, transpose, levels_[0].bands[i],
                  levels_[0].bands[i + 1]);
    }));
  }
  tasks.wait();
  tasks.get();
}

void Pyramid::BlurXThread(float radius, Pyramid* transpose, int sy, int ey) {
  int x;
  int y;
  int i;
  int o;
  float* line0 =
      levels_[0].data.get() + static_cast<ptrdiff_t>(sy) * levels_[0].pitch;
  float* line1 = line0 + levels_[0].pitch;
  float* line2 = line1 + levels_[0].pitch;
  float* line3 = line2 + levels_[0].pitch;
  float* out = transpose->levels_[0].data.get() + sy;
  simde__m128 temp1;
  simde__m128 temp2;

  int iradius = (int)std::floor(radius);
  simde__m128 irp1 = simde_mm_set_ps1((float)(iradius + 1));
  simde__m128 mul = simde_mm_set_ps1(radius - iradius);
  simde__m128 acc;

  int left;
  int right;

  auto blur_sse_get = [line3, line2, line1, line0](int x) {
    return simde_mm_set_ps(line3[x], line2[x], line1[x], line0[x]);
  };

  // +3 is probably not necessary because all bands are mod 4
  int fours = (ey - sy + 3) >> 2;

  if (iradius < levels_[0].width >> 1) {
    for (y = 0; y < fours; ++y) {
      acc = simde_mm_setzero_ps();
      left = 0;

      temp1 = blur_sse_get(left);
      left++;

      acc = simde_mm_mul_ps(temp1, irp1);
      for (right = 1; right < iradius + 1;) {
        temp2 = blur_sse_get(right);
        right++;
        acc = simde_mm_add_ps(acc, temp2);
      }

      x = 0;
      o = 0;
      right = iradius + 1;

      for (i = 0; i <= iradius; ++i) {
        temp2 = blur_sse_get(right);
        right++;
        simde_mm_store_ps(
            &out[o],
            simde_mm_add_ps(
                acc, simde_mm_mul_ps(simde_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        acc = simde_mm_add_ps(simde_mm_sub_ps(temp2, temp1), acc);
        ++x;
      }

      while (right < levels_[0].width) {
        temp2 = blur_sse_get(right);
        right++;
        simde_mm_store_ps(
            &out[o],
            simde_mm_add_ps(
                acc, simde_mm_mul_ps(simde_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        temp1 = blur_sse_get(left);
        left++;
        acc = simde_mm_add_ps(simde_mm_sub_ps(temp2, temp1), acc);
        ++x;
      }

      while (x < levels_[0].width) {
        simde_mm_store_ps(
            &out[o],
            simde_mm_add_ps(
                acc, simde_mm_mul_ps(simde_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        temp1 = blur_sse_get(left);
        left++;
        acc = simde_mm_add_ps(simde_mm_sub_ps(temp2, temp1), acc);
        ++x;
      }

      line0 += levels_[0].pitch << 2;
      line1 += levels_[0].pitch << 2;
      line2 += levels_[0].pitch << 2;
      line3 += levels_[0].pitch << 2;
      out += 4;
    }
  } else {
    // if radius is wider than image
    for (y = 0; y < fours; ++y) {
      acc = simde_mm_setzero_ps();

      temp1 = blur_sse_get(0);
      acc = simde_mm_mul_ps(temp1, irp1);
      right = 1;
      for (x = 1; x < iradius + 1; ++x) {
        if (right < levels_[0].width) {
          temp2 = blur_sse_get(right);
          ++right;
        }
        acc = simde_mm_add_ps(acc, temp2);
      }

      x = 0;
      o = 0;
      left = -iradius;

      for (x = 0; x < levels_[0].width; ++x) {
        if (right < levels_[0].width) {
          temp2 = blur_sse_get(right);
          ++right;
        }
        simde_mm_store_ps(
            &out[o],
            simde_mm_add_ps(
                acc, simde_mm_mul_ps(simde_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        if (left > 0) {
          temp1 = blur_sse_get(left);
        }
        ++left;
        acc = simde_mm_add_ps(simde_mm_sub_ps(temp2, temp1), acc);
      }

      line0 += levels_[0].pitch << 2;
      line1 += levels_[0].pitch << 2;
      line2 += levels_[0].pitch << 2;
      line3 += levels_[0].pitch << 2;
      out += 4;
    }
  }
}

/***********************************************************************
 * PNG debug
 ***********************************************************************/
void Pyramid::Png(const char* filename) {
  int width = levels_[0].pitch;
  int height =
      levels_[0].height + (levels_.size() > 1 ? 1 + levels_[1].height : 0);

  auto temp = std::vector<uint8_t>(static_cast<std::size_t>(width) * height, 0);

  int px = 0;
  int py = 0;

  for (int l = 0; l < (int)levels_.size(); ++l) {
    auto* data = (float*)levels_[l].data.get();
    uint8_t* line =
        temp.data() + static_cast<ptrdiff_t>(py) * levels_[0].pitch + px;
    for (int y = 0; y < levels_[l].height; ++y) {
      for (int x = 0; x < levels_[l].pitch; ++x) {
        int f = (int)floor(data[x] + 0.5);
        line[x] = std::max(0, std::min(255, f));
      }
      line += levels_[0].pitch;
      data += levels_[l].pitch;
    }
    if ((l & 1) != 0) {
      px += levels_[l].pitch + 1;
    } else {
      py += levels_[l].height + 1;
    }
  }

  io::png::Pnger::Quick(filename, temp.data(), width, height, width,
                        io::png::ColorType::GRAY);
}

/***********************************************************************
************************************************************************
* Output
************************************************************************
***********************************************************************/
namespace {
/***********************************************************************
 * 8-bit planar
 ***********************************************************************/
template <typename F>
void OutPlanar8(void* dst_vp, F loader, int pitch, const Pyramid::Level& level,
                int band, bool chroma) {
  int x;
  int y;
  auto* dst_p = (uint8_t*)dst_vp;
  uint8_t black = chroma ? 0x80 : 0x00;

  simde__m128 zeroes = simde_mm_setzero_ps();
  simde__m128 maxes = simde_mm_set_ps1(255.0f);

  simde__m128i shuffle1 =
      simde_mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x0c080400);
  simde__m128i shuffle2 =
      simde_mm_set_epi32(0x80808080, 0x80808080, 0x0c080400, 0x80808080);
  simde__m128i shuffle3 =
      simde_mm_set_epi32(0x80808080, 0x0c080400, 0x80808080, 0x80808080);
  simde__m128i shuffle4 =
      simde_mm_set_epi32(0x0c080400, 0x80808080, 0x80808080, 0x80808080);
  simde__m128i four_shuffle =
      simde_mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x0c080400);
  simde__m128i pixels;
  simde__m128* p_p;
  auto* p_pt = (simde__m128*)level.data.get();

  simde__m128i* dst_pp_m;
  int* dst_pp_i;
  uint8_t* dst_pp_b;

  int m128_pitch = level.pitch >> 2;

  int sixteens = level.width >> 4;
  int fours = (level.width >> 2) - (sixteens << 2);
  int singles = level.width & 3;

  int sy = level.bands[band];
  int ey = level.bands[band + 1];

  if (level.id > 0) {
    if (sy == 0) {
      sy++;
    }
    if (ey == level.height) {
      ey--;
    }
  }

  dst_p += (std::size_t)(sy - (level.id ? 1 : 0)) * pitch;

  p_pt += (std::size_t)m128_pitch * sy;

  simde__m128 dither_add;

  auto load = [&]() -> simde__m128 {
    return loader((float*)p_p++, dither_add, zeroes, maxes);
  };

  for (y = sy; y < ey; y++) {
    switch (y & 3) {
      case 0:
        dither_add = simde_mm_set_ps(0.4999f, 0.0f, 0.375f, -0.125f);
        break;
      case 1:
        dither_add = simde_mm_set_ps(-0.25f, 0.25f, -0.375f, 0.125f);
        break;
      case 2:
        dither_add = simde_mm_set_ps(0.3125f, -0.1875f, 0.4375f, -0.0625f);
        break;
      case 3:
        dither_add = simde_mm_set_ps(-0.4375f, 0.0625f, -0.3125f, 0.1875f);
        break;
    }

    dst_pp_m = (simde__m128i*)dst_p;
    p_p = p_pt;
    for (x = 0; x < sixteens; ++x) {
      pixels = simde_mm_shuffle_epi8(simde_mm_cvtps_epi32(load()), shuffle1);
      pixels = simde_mm_or_si128(
          pixels,
          simde_mm_shuffle_epi8(simde_mm_cvtps_epi32(load()), shuffle2));
      pixels = simde_mm_or_si128(
          pixels,
          simde_mm_shuffle_epi8(simde_mm_cvtps_epi32(load()), shuffle3));
      pixels = simde_mm_or_si128(
          pixels,
          simde_mm_shuffle_epi8(simde_mm_cvtps_epi32(load()), shuffle4));

      simde_mm_storeu_si128(dst_pp_m++, pixels);
    }

    dst_pp_i = (int*)dst_pp_m;
    for (x = 0; x < fours; ++x) {
      *dst_pp_i++ = simde_mm_cvtsi128_si32(
          simde_mm_shuffle_epi8(simde_mm_cvtps_epi32(load()), four_shuffle));
    }

    if (singles) {
      dst_pp_b = (uint8_t*)dst_pp_i;
      simde__m128i a;
      a = simde_mm_cvtps_epi32(load());
      *dst_pp_b++ = (uint8_t)simde_mm_extract_epi8(a, 0);
      if (singles > 1) {
        *dst_pp_b++ = (uint8_t)simde_mm_extract_epi8(a, 4);
        if (singles == 3) {
          *dst_pp_b++ = (uint8_t)simde_mm_extract_epi8(a, 8);
        }
      }
    }

    if (level.id > 0) {
      memcpy(dst_p, dst_p + 1, level.width - 2);
      dst_p[level.width - 1] = (uint8_t)black;
      dst_p[level.width - 2] = (uint8_t)black;
    }

    p_pt += m128_pitch;
    dst_p += pitch;
  }
};

/***********************************************************************
 * 16-bit planar
 ***********************************************************************/
template <typename F>
void OutPlanar16(void* dst_vp, F loader, int pitch, const Pyramid::Level& level,
                 int band, bool chroma) {
  int x;
  int y;
  auto* dst_p = (uint16_t*)dst_vp;
  uint16_t black = chroma ? 0x8000 : 0x0000;

  simde__m128 zeroes = simde_mm_setzero_ps();
  simde__m128 maxes = simde_mm_set_ps1(65535.0f);

  simde__m128i shuffle1 =
      simde_mm_set_epi32(0x80808080, 0x80808080, 0x0d0c0908, 0x05040100);
  simde__m128i shuffle2 =
      simde_mm_set_epi32(0x0d0c0908, 0x05040100, 0x80808080, 0x80808080);
  simde__m128i pixels;
  simde__m128* p_p;
  auto* p_pt = (simde__m128*)level.data.get();

  simde__m128i* dst_pp_m;
  uint16_t* dst_pp_w;

  int m128_pitch = level.pitch >> 2;

  int eights = level.width >> 3;
  int four = level.width & 4;
  int singles = level.width & 3;

  int sy = level.bands[band];
  int ey = level.bands[band + 1];

  if (level.id > 0) {
    if (sy == 0) {
      sy++;
    }
    if (ey == level.height) {
      ey--;
    }
  }

  dst_p += (static_cast<ptrdiff_t>(sy - (level.id ? 1 : 0))) * pitch;

  p_pt += static_cast<ptrdiff_t>(m128_pitch) * sy;

  simde__m128 dither_add;

  auto load = [&]() -> simde__m128 {
    return loader((float*)p_p++, dither_add, zeroes, maxes);
  };

  for (y = sy; y < ey; y++) {
    switch (y & 3) {
      case 0:
        dither_add = simde_mm_set_ps(0.4999f, 0.0f, 0.375f, -0.125f);
        break;
      case 1:
        dither_add = simde_mm_set_ps(-0.25f, 0.25f, -0.375f, 0.125f);
        break;
      case 2:
        dither_add = simde_mm_set_ps(0.3125f, -0.1875f, 0.4375f, -0.0625f);
        break;
      case 3:
        dither_add = simde_mm_set_ps(-0.4375f, 0.0625f, -0.3125f, 0.1875f);
        break;
    }

    dst_pp_m = (simde__m128i*)dst_p;
    p_p = p_pt;
    for (x = 0; x < eights; ++x) {
      pixels = simde_mm_shuffle_epi8(simde_mm_cvtps_epi32(load()), shuffle1);
      pixels = simde_mm_or_si128(
          pixels,
          simde_mm_shuffle_epi8(simde_mm_cvtps_epi32(load()), shuffle2));

      simde_mm_storeu_si128(dst_pp_m++, pixels);
    }

    simde__m128i a;

    dst_pp_w = (uint16_t*)dst_pp_m;
    if (four) {
      a = simde_mm_cvtps_epi32(load());
      *dst_pp_w++ = (uint16_t)simde_mm_extract_epi16(a, 0);
      *dst_pp_w++ = (uint16_t)simde_mm_extract_epi16(a, 2);
      *dst_pp_w++ = (uint16_t)simde_mm_extract_epi16(a, 4);
      *dst_pp_w++ = (uint16_t)simde_mm_extract_epi16(a, 6);
    }

    if (singles) {
      a = simde_mm_cvtps_epi32(load());
      *dst_pp_w++ = (uint16_t)simde_mm_extract_epi16(a, 0);
      if (singles > 1) {
        *dst_pp_w++ = (uint16_t)simde_mm_extract_epi16(a, 2);
        if (singles == 3) {
          *dst_pp_w++ = (uint16_t)simde_mm_extract_epi16(a, 4);
        }
      }
    }

    if (level.id > 0) {
      memcpy(dst_p, dst_p + 1, (level.width - 2) << 1);
      dst_p[level.width - 1] = black;
      dst_p[level.width - 2] = black;
    }

    p_pt += m128_pitch;
    dst_p += pitch;
  }
}

/***********************************************************************
 * 32-bit planar
 ***********************************************************************/

template <typename F>
void OutPlanar32(void* dst_vp, F loader, int pitch, const Pyramid::Level& level,
                 int band, bool chroma) {
  int x;
  int y;
  auto* dst_p = (simde__m128*)dst_vp;

  simde__m128 zeroes;
  simde__m128 maxes;
  simde__m128 dither_add = simde_mm_set_ps1(0.0f);

  pitch >>= 2;  // number of floats to number of simde__m128s

  if (chroma) {
    zeroes = simde_mm_set_ps1(-0.5f);
    maxes = simde_mm_set_ps1(0.5f);
  } else {
    zeroes = simde_mm_set_ps1(0.0f);
    maxes = simde_mm_set_ps1(1.0f);
  }

  auto* p_p = (simde__m128*)level.data.get();

  auto load = [&]() -> simde__m128 {
    return loader((float*)&p_p[x], dither_add, zeroes, maxes);
  };

  float* dst_pp_f;

  int m128_pitch = level.pitch >> 2;

  int fours = (level.width + 3) >> 2;
  int wipes = (fours << 2) - level.width;

  int sy = level.bands[band];
  int ey = level.bands[band + 1];

  if (level.id > 0) {
    if (sy == 0) {
      sy++;
    }
    if (ey == level.height) {
      ey--;
    }
  }

  dst_p += (static_cast<ptrdiff_t>(sy - (level.id ? 1 : 0))) * pitch;

  p_p += static_cast<ptrdiff_t>(m128_pitch) * sy;

  for (y = sy; y < ey; y++) {
    for (x = 0; x < fours; x++) {
      dst_p[x] = load();
    }

    if (wipes) {
      dst_pp_f = (float*)&dst_p[x];
      int w = wipes;
      while (w) {
        *--dst_pp_f = 0;
        w--;
      }
    }

    if (level.id > 0) {
      memcpy(dst_p, ((float*)dst_p) + 1, (level.width - 2) << 2);
      ((float*)dst_p)[level.width - 1] = 0.0f;
      ((float*)dst_p)[level.width - 2] = 0.0f;
    }

    p_p += m128_pitch;
    dst_p += pitch;
  }
}

/***********************************************************************
 * Interleaved
 ***********************************************************************/
template <typename T, typename F>
void OutInterleaved(T dst_p, F loader, int pitch, const Pyramid::Level& level,
                    int band, bool chroma, int step, int offset) {
  int x;
  int y;

  simde__m128 zeroes = simde_mm_setzero_ps();
  simde__m128 maxes = simde_mm_set_ps1(sizeof(*dst_p) == 1 ? 255.0f : 65535.0f);

  int sy = level.bands[band];
  int ey = level.bands[band + 1];

  if (level.id > 0) {
    if (sy == 0) {
      sy++;
    }
    if (ey == level.height) {
      ey--;
    }
  }

  dst_p += (sy - (level.id ? 1 : 0)) * pitch + offset;

  int m128_pitch = level.pitch >> 2;
  simde__m128* p_p =
      (simde__m128*)level.data.get() + static_cast<ptrdiff_t>(m128_pitch) * sy;
  T dst_pp;

  simde__m128i a;
  int fours = level.width >> 2;
  int singles = level.width & 3;

  if (level.id > 0) {
    singles--;
    if (singles < 0) {
      singles = 3;
      fours--;
    }
  }

  simde__m128 dither_add;

  auto load = [&]() -> simde__m128 {
    return loader((float*)&p_p[x], dither_add, zeroes, maxes);
  };

  // loop
  for (y = sy; y < ey; y++) {
    switch (y & 3) {
      case 0:
        dither_add = simde_mm_set_ps(0.4999f, 0.0f, 0.375f, -0.125f);
        break;
      case 1:
        dither_add = simde_mm_set_ps(-0.25f, 0.25f, -0.375f, 0.125f);
        break;
      case 2:
        dither_add = simde_mm_set_ps(0.3125f, -0.1875f, 0.4375f, -0.0625f);
        break;
      case 3:
        dither_add = simde_mm_set_ps(-0.4375f, 0.0625f, -0.3125f, 0.1875f);
        break;
    }

    dst_pp = dst_p;
    x = 0;
    if (level.id > 0) {
      a = simde_mm_cvtps_epi32(load());
      *dst_pp = simde_mm_extract_epi16(a, 2);
      dst_pp += step;
      *dst_pp = simde_mm_extract_epi16(a, 4);
      dst_pp += step;
      *dst_pp = simde_mm_extract_epi16(a, 6);
      dst_pp += step;
      x++;
    }

    for (; x < fours; x++) {
      a = simde_mm_cvtps_epi32(load());
      *dst_pp = simde_mm_extract_epi16(a, 0);
      dst_pp += step;
      *dst_pp = simde_mm_extract_epi16(a, 2);
      dst_pp += step;
      *dst_pp = simde_mm_extract_epi16(a, 4);
      dst_pp += step;
      *dst_pp = simde_mm_extract_epi16(a, 6);
      dst_pp += step;
    }

    if (singles) {
      a = simde_mm_cvtps_epi32(load());
      *dst_pp = simde_mm_extract_epi8(a, 0);
      if (singles > 1) {
        dst_pp += step;
        *dst_pp = simde_mm_extract_epi16(a, 2);
        if (singles == 3) {
          dst_pp += step;
          *dst_pp = simde_mm_extract_epi16(a, 2);
        }
      }
    }

    p_p += m128_pitch;
    dst_p += pitch;
  }
}

/***********************************************************************
 * Loaders
 ***********************************************************************/

enum class Transform {
  kNone,
  kGamma,
  kDither,
  kClamp,
  kDitherGamma,
  kClampGamma,
  kClampDither,
  kClampDitherGamma,
};

template <Transform mode>
struct Loader {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    return simde_mm_load_ps(src_p);
  }
};

template <>
struct Loader<Transform::kGamma> {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    return simde_mm_sqrt_ps(simde_mm_load_ps(src_p));
  }
};

template <>
struct Loader<Transform::kDither> {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    return simde_mm_add_ps(simde_mm_load_ps(src_p), dither_add);
  }
};

template <>
struct Loader<Transform::kClamp> {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    return simde_mm_min_ps(simde_mm_max_ps(simde_mm_load_ps(src_p), zeroes),
                           maxes);
  }
};

template <>
struct Loader<Transform::kDitherGamma> {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    return simde_mm_add_ps(simde_mm_sqrt_ps(simde_mm_load_ps(src_p)),
                           dither_add);
  }
};

template <>
struct Loader<Transform::kClampGamma> {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    // Originally:
    // simde_mm_min_ps(simde_mm_sqrt_ps(simde_mm_max_ps(simde_mm_load_ps(p),
    // z)), m)
    return simde_mm_min_ps(
        simde_mm_max_ps(simde_mm_sqrt_ps(simde_mm_load_ps(src_p)), zeroes),
        maxes);
  }
};

template <>
struct Loader<Transform::kClampDither> {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    // Originally:
    // simde_mm_min_ps(simde_mm_add_ps(simde_mm_max_ps(simde_mm_load_ps(p), z),
    // d), m)
    return simde_mm_min_ps(
        simde_mm_max_ps(simde_mm_add_ps(simde_mm_load_ps(src_p), dither_add),
                        zeroes),
        maxes);
  }
};

template <>
struct Loader<Transform::kClampDitherGamma> {
  simde__m128 operator()(float* src_p, simde__m128 dither_add,
                         simde__m128 zeroes, simde__m128 maxes) {
    // Originally:
    // simde_mm_min_ps(simde_mm_add_ps(simde_mm_sqrt_ps(simde_mm_max_ps(simde_mm_load_ps(p),
    // z)), d), m)
    return simde_mm_min_ps(
        simde_mm_max_ps(
            simde_mm_add_ps(simde_mm_sqrt_ps(simde_mm_load_ps(src_p)),
                            dither_add),
            zeroes),
        maxes);
  }
};

}  // namespace

/***********************************************************************
 * Out
 ***********************************************************************/
template <typename T>
void Pyramid::Out(T dst_p, int pitch, bool gamma, bool dither, bool clamp,
                  int level, int step, int offset, bool chroma) {
  int bytes = sizeof(*dst_p);
  int eb = (int)(levels_[level].bands.size() - 1);

  using Type =
      typename std::conditional<sizeof(*dst_p) == 1, uint8_t*, uint16_t*>::type;
  // used to avoid generating a float* version of OutInterleaved which would
  // cause warnings

  int s = (gamma ? 1 : 0) | (dither && bytes != 4 ? 2 : 0) | (clamp ? 4 : 0);

  auto tasks = mt::MultiFuture<void>{};
  if (step) {  // interleaved
    for (int band = 0; band < eb; ++band) {
      switch (s) {
        case 0:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kNone>{}, pitch,
                           levels_[level], band, chroma, step, offset);
          }));
          break;
        case 1:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kGamma>{}, pitch,
                           levels_[level], band, chroma, step, offset);
          }));
          break;
        case 2:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kDither>{}, pitch,
                           levels_[level], band, chroma, step, offset);
          }));
          break;
        case 3:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kDitherGamma>{},
                           pitch, levels_[level], band, chroma, step, offset);
          }));
          break;
        case 4:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kClamp>{}, pitch,
                           levels_[level], band, chroma, step, offset);
          }));
          break;
        case 5:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kClampGamma>{}, pitch,
                           levels_[level], band, chroma, step, offset);
          }));
          break;
        case 6:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kClampDither>{},
                           pitch, levels_[level], band, chroma, step, offset);
          }));
          break;
        case 7:
          tasks.push_back(threadpool_->Queue([=, this] {
            OutInterleaved((Type)dst_p, Loader<Transform::kClampDitherGamma>{},
                           pitch, levels_[level], band, chroma, step, offset);
          }));
          break;
      }
    }
  } else {  // planar
    switch (bytes) {
      case 1: {
        for (int band = 0; band < eb; ++band) {
          switch (s) {
            case 0:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kNone>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
            case 1:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kGamma>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
            case 2:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kDither>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
            case 3:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kDitherGamma>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
            case 4:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kClamp>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
            case 5:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kClampGamma>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
            case 6:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kClampDither>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
            case 7:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar8(dst_p, Loader<Transform::kClampDitherGamma>{}, pitch,
                           levels_[level], band, chroma);
              }));
              break;
          }
        }
      } break;
      case 2: {
        for (int band = 0; band < eb; ++band) {
          switch (s) {
            case 0:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kNone>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 1:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kGamma>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 2:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kDither>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 3:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kDitherGamma>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 4:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kClamp>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 5:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kClampGamma>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 6:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kClampDither>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 7:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar16(dst_p, Loader<Transform::kClampDitherGamma>{},
                            pitch, levels_[level], band, chroma);
              }));
              break;
          }
        }
      } break;
      case 4: {
        for (int band = 0; band < eb; ++band) {
          switch (s) {
            case 0:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar32(dst_p, Loader<Transform::kNone>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 1:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar32(dst_p, Loader<Transform::kGamma>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 4:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar32(dst_p, Loader<Transform::kClamp>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
            case 5:
              tasks.push_back(threadpool_->Queue([=, this] {
                OutPlanar32(dst_p, Loader<Transform::kClampGamma>{}, pitch,
                            levels_[level], band, chroma);
              }));
              break;
          }
        }
      } break;
    }
  }

  tasks.wait();
  tasks.get();
}

template void Pyramid::Out(uint16_t* dst_p, int pitch, bool gamma, bool dither,
                           bool clamp, int level, int step, int offset,
                           bool chroma);

template void Pyramid::Out(uint8_t* dst_p, int pitch, bool gamma, bool dither,
                           bool clamp, int level, int step, int offset,
                           bool chroma);

}  // namespace multiblend
