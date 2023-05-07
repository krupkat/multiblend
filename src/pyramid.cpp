#include "src/pyramid.h"

#ifdef PYR_DENOISE
#define USE_SSE2
#include "sse_mathfun.h"
#endif

#include <math.h>

#include "src/linux_overrides.h"
#include "src/pnger.h"
#include "src/threadpool.h"

namespace multiblend {

/***********************************************************************
 * Constructor/destructor
 ***********************************************************************/
Pyramid::Pyramid(int width, int height, int _levels, Pyramid* share)
    : Pyramid(width, height, _levels, 0, 0, false, share) {}

Pyramid::Pyramid(int width, int height, int _levels, int x, int y,
                 bool _no_alloc, Pyramid* share)
    : shared_(share != NULL), no_alloc_(_no_alloc) {
  float* data;

  int n_levels = DefaultNumLevels(width, height);
  if (_levels) {
    if (_levels < 0)
      n_levels -= _levels;
    else
      n_levels = _levels;
  }

  int b, req_alignment;

  b = 0;
  req_alignment = 2;

  for (int n = 0; n < n_levels; ++n) {
    bool x_shift = (x - b) & (req_alignment - 1);
    bool y_shift = (y - b) & (req_alignment - 1);
    int pitch = (width + x_shift + 7) & ~7;
    size_t bytes =
        (size_t)pitch * ((height + y_shift + 3) & ~3) * sizeof(float);
    total_bytes_ += bytes;

    if (shared_) {
      data = share->levels_[levels_.size()].data;
    } else {
      if (!no_alloc_) {
        data = (float*)_aligned_malloc(bytes, 16);  // was 32
        if (!data) {
          for (int j = 0; j < n; ++j) _aligned_free(levels_[j].data);
          throw(bytes);
        }
      } else {
        data = NULL;
      }
    }

    levels_.push_back({width, height, pitch, pitch >> 2, bytes, data, x, y,
                       x_shift, y_shift, n ? levels_[n - 1].x_shift : false,
                       n ? levels_[n - 1].m128_pitch : 0});

    x -= (x_shift << n);
    y -= (y_shift << n);

    x -= req_alignment;
    y -= req_alignment;

    b -= req_alignment;

    req_alignment <<= 1;

    width = (width + x_shift + 6) >> 1;
    height = (height + y_shift + 6) >> 1;
  }

  threadpool_ = mt::Threadpool::GetInstance();

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

  lines_ = new __m128*[threadpool_->GetNThreads()];

  for (int i = 0; i < threadpool_->GetNThreads(); ++i) {
    lines_[i] = (__m128*)_aligned_malloc(levels_[0].pitch * sizeof(float),
                                         16);  // was 32
  }
}

Pyramid::~Pyramid() {
  for (int i = 0; i < threadpool_->GetNThreads(); ++i) {
    _aligned_free(lines_[i]);
  }
  delete lines_;

  if (!shared_ && !no_alloc_) {
    for (auto it = levels_.begin(); it < levels_.end(); ++it) {
      _aligned_free(it->data);
    }
  }
  levels_.clear();

  free(lut_);
}

/***********************************************************************
 * copiers
 ***********************************************************************/
void Pyramid::set_lut(int bits, bool gamma) {
  if (lut_bits_ < bits || lut_gamma_ != gamma || !lut_) {
    free(lut_);
    lut_ = (float*)malloc((1 << bits) << 2);
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

    for (int t = 0; t < (int)levels_[0].bands.size() - 1; ++t) {
      switch (bits) {
        case 8:
          threadpool_->Queue([=, this] {
            CopyInterleavedThread_8bit(src_p, step, pitch, levels_[0].bands[t],
                                       levels_[0].bands[t + 1]);
          });
          break;
        case 16:
          threadpool_->Queue([=, this] {
            CopyInterleavedThread_16bit((uint16_t*)src_p, step, pitch,
                                        levels_[0].bands[t],
                                        levels_[0].bands[t + 1]);
          });
          break;
        case 32:
          break;
      }
    }
    threadpool_->Wait();
  } else {
    for (int t = 0; t < (int)levels_[0].bands.size() - 1; ++t) {
      // planar (only slight improvement with MT, but increases CPU usage)
      switch (bits) {
        case 8:
          threadpool_->Queue([=, this] {
            CopyPlanarThread_8bit(src_p, pitch, gamma, levels_[0].bands[t],
                                  levels_[0].bands[t + 1]);
          });
          break;
        case 10:
        case 12:
        case 14:
        case 16:
          threadpool_->Queue([=, this] {
            CopyPlanarThread_16bit((uint16_t*)src_p, pitch, gamma,
                                   levels_[0].bands[t],
                                   levels_[0].bands[t + 1]);
          });
          break;
        case 32:
          threadpool_->Queue([=, this] {
            CopyPlanarThread_32bit((__m128*)src_p, pitch, gamma,
                                   levels_[0].bands[t],
                                   levels_[0].bands[t + 1]);
          });
          break;
      }
    }

    threadpool_->Wait();
  }
}

void Pyramid::CopyInterleavedThread_8bit(uint8_t* src_p, int step, int pitch,
                                         int sy, int ey) {
  int x, y;
  uint8_t* src_pp;

  src_p += sy * pitch;

  float* p_p = levels_[0].data;
  p_p += sy * levels_[0].pitch;

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
  int x, y;
  uint16_t* src_pp;

  src_p += sy * pitch;

  float* p_p = levels_[0].data;
  p_p += sy * levels_[0].pitch;

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
  int x, y;
  __m128i pixels;
  __m128 fpixels;
  __m128i shuffle =
      _mm_set_epi32(0x80808003, 0x80808002, 0x80808001, 0x80808000);
  __m128* rp_p = (__m128*)levels_[0].data;
  __m128* p_p;
  __m128i* src_pp_m;
  int* src_pp_i;
  uint8_t* src_pp_b;

  int sixteens = levels_[0].width >> 4;
  int fours = (levels_[0].width & 0xf) >> 2;
  int ones = (levels_[0].width & 3);
  int extras = levels_[0].pitch - levels_[0].width;

  src_p += sy * pitch;
  rp_p += sy * levels_[0].m128_pitch;

  int g;

  for (y = sy; y < ey; ++y) {
    src_pp_m = (__m128i*)src_p;
    p_p = rp_p;
    if (gamma) {
      for (x = 0; x < sixteens; ++x) {
        pixels = _mm_loadu_si128(src_pp_m++);
        fpixels = _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle));
        _mm_store_ps((float*)p_p++, _mm_mul_ps(fpixels, fpixels));

        pixels = _mm_srli_si128(pixels, 4);
        fpixels = _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle));
        _mm_store_ps((float*)p_p++, _mm_mul_ps(fpixels, fpixels));

        pixels = _mm_srli_si128(pixels, 4);
        fpixels = _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle));
        _mm_store_ps((float*)p_p++, _mm_mul_ps(fpixels, fpixels));

        pixels = _mm_srli_si128(pixels, 4);
        fpixels = _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle));
        _mm_store_ps((float*)p_p++, _mm_mul_ps(fpixels, fpixels));
      }
    } else {
      for (x = 0; x < sixteens; ++x) {
        pixels = _mm_loadu_si128(src_pp_m++);
        _mm_store_ps((float*)p_p++,
                     _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle)));

        pixels = _mm_srli_si128(pixels, 4);
        _mm_store_ps((float*)p_p++,
                     _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle)));

        pixels = _mm_srli_si128(pixels, 4);
        _mm_store_ps((float*)p_p++,
                     _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle)));

        pixels = _mm_srli_si128(pixels, 4);
        _mm_store_ps((float*)p_p++,
                     _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle)));
      }
    }

    src_pp_i = (int*)src_pp_m;

    for (x = 0; x < fours; ++x) {
      if (gamma) {
        fpixels = _mm_cvtepi32_ps(
            _mm_shuffle_epi8(_mm_cvtsi32_si128(*src_pp_i++), shuffle));
        _mm_store_ps((float*)p_p++, _mm_mul_ps(fpixels, fpixels));
      } else {
        _mm_store_ps((float*)p_p++,
                     _mm_cvtepi32_ps(_mm_shuffle_epi8(
                         _mm_cvtsi32_si128(*src_pp_i++), shuffle)));
      }
    }

    src_pp_b = (uint8_t*)src_pp_i;
    float* fp_p = (float*)p_p;

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
  int x, y;

  __m128i pixels;
  __m128 fpixels;
  __m128i shuffle1 =
      _mm_set_epi32(0x80800706, 0x80800504, 0x80800302, 0x80800100);
  __m128i shuffle2 =
      _mm_set_epi32(0x80800f0e, 0x80800d0c, 0x80800b0a, 0x80800908);
  __m128* rp_p = (__m128*)levels_[0].data;
  __m128* p_p;
  __m128i* src_pp_m;
  uint16_t* src_pp_w;

  int eights = levels_[0].width >> 3;
  int ones = (levels_[0].width & 7);
  int extras = levels_[0].pitch - levels_[0].width;

  src_p += sy * pitch;
  rp_p += sy * levels_[0].m128_pitch;

  int g;

  for (y = sy; y < ey; ++y) {
    src_pp_m = (__m128i*)src_p;
    p_p = rp_p;
    if (gamma) {
      for (x = 0; x < eights; ++x) {
        pixels = _mm_loadu_si128(src_pp_m++);
        fpixels = _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle1));
        _mm_store_ps((float*)p_p++, _mm_mul_ps(fpixels, fpixels));

        fpixels = _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle2));
        _mm_store_ps((float*)p_p++, _mm_mul_ps(fpixels, fpixels));
      }
    } else {
      for (x = 0; x < eights; ++x) {
        pixels = _mm_loadu_si128(src_pp_m++);
        _mm_store_ps((float*)p_p++,
                     _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle1)));

        _mm_store_ps((float*)p_p++,
                     _mm_cvtepi32_ps(_mm_shuffle_epi8(pixels, shuffle2)));
      }
    }

    src_pp_w = (uint16_t*)src_pp_m;
    float* fp_p = (float*)p_p;

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

void Pyramid::CopyPlanarThread_32bit(__m128* src_p, int pitch, bool gamma,
                                     int sy, int ey) {
  int x, y;

  pitch >>= 2;

  __m128* rp_p = (__m128*)levels_[0].data;

  src_p += sy * pitch;
  rp_p += sy * levels_[0].m128_pitch;

  if (gamma) {
    int fours = (levels_[0].width + 3) >> 2;

    for (y = sy; y < ey; ++y) {
      for (x = 0; x < fours; ++x) {
        __m128 pixels = _mm_load_ps((float*)&src_p[x]);
        _mm_store_ps((float*)&rp_p[x], _mm_mul_ps(pixels, pixels));
      }
      float copy = ((float*)rp_p)[levels_[0].width - 1];
      for (x = levels_[0].width; x < levels_[0].pitch; ++x)
        ((float*)rp_p)[x] = copy;
      src_p += pitch;
      rp_p += levels_[0].m128_pitch;
    }
  } else {
    int copy_bytes = levels_[0].width << 2;

    for (y = sy; y < ey; ++y) {
      memcpy(rp_p, src_p, copy_bytes);
      float copy = ((float*)rp_p)[levels_[0].width - 1];
      for (x = levels_[0].width; x < levels_[0].pitch; ++x)
        ((float*)rp_p)[x] = copy;
      src_p += pitch;
      rp_p += levels_[0].m128_pitch;
    }
  }
}

/***********************************************************************
 * subsample
 ***********************************************************************/
void Pyramid::Subsample(int sub_w, int sub_h, Pyramid* source) {
  int x, y;
  int p = 0;
  __m128* in = (__m128*)source->levels_[0].data;
  __m128* out = (__m128*)levels_[0].data;
  __m128** temp_lines = source->lines_;
  __m128* line = temp_lines[0];
  int m128_pitch_in = source->levels_[0].m128_pitch;
  int m128_pitch_out = levels_[0].m128_pitch;
  int mid_pitch = sub_w == 2 ? ((m128_pitch_in >> 1) + 1) & ~1 : m128_pitch_in;
  __m128 three = _mm_set_ps1(3);
  __m128 four = _mm_set_ps1(4);
  __m128 mul = _mm_set_ps1((sub_h ? 1.0f / 8 : 1.0f) *
                           (sub_w == 2 ? 1.0f / 64 : 1.0f / 8));

  for (y = 0; y < levels_[0].height; ++y) {
    if (sub_h) {
      if (y == 0) {
        for (x = 0; x < m128_pitch_in; ++x) {
          _mm_store_ps(
              (float*)&temp_lines[0][x],
              _mm_add_ps(
                  _mm_mul_ps(_mm_load_ps((float*)&in[x]), four),
                  _mm_add_ps(
                      _mm_mul_ps(_mm_load_ps((float*)&in[x + m128_pitch_in]),
                                 three),
                      _mm_load_ps((float*)&in[x + (m128_pitch_in << 1)]))));
        }
      } else if (y == levels_[0].height - 1) {
        for (x = 0; x < m128_pitch_in; ++x) {
          _mm_store_ps(
              (float*)&temp_lines[0][x],
              _mm_add_ps(
                  _mm_mul_ps(_mm_load_ps((float*)&in[x + m128_pitch_in]), four),
                  _mm_add_ps(_mm_mul_ps(_mm_load_ps((float*)&in[x]), three),
                             _mm_load_ps((float*)&in[x - m128_pitch_in]))));
        }
      } else {
        for (x = 0; x < m128_pitch_in; ++x) {
          _mm_store_ps(
              (float*)&temp_lines[0][x],
              _mm_add_ps(
                  _mm_add_ps(
                      _mm_load_ps((float*)&in[x - m128_pitch_in]),
                      _mm_load_ps((float*)&in[x + (m128_pitch_in << 1)])),
                  _mm_mul_ps(
                      _mm_add_ps(_mm_load_ps((float*)&in[x]),
                                 _mm_load_ps((float*)&in[x + m128_pitch_in])),
                      three)));
        }
      }
      in += m128_pitch_in << 1;
    } else {
      line = (__m128*)in;
      in += m128_pitch_in;
    }
    switch (sub_w) {
      case 2:
        Subsample_Squeeze(line, line, m128_pitch_in, mid_pitch, NULL);
      case 1:
        Subsample_Squeeze(line, out, mid_pitch, m128_pitch_out, &mul);
    }
    out += m128_pitch_out;
  }
}

void Pyramid::Subsample_Squeeze(__m128* in, __m128* Out, int m128_pitch_in,
                                int m128_pitch_out, __m128* mul) {
  int read = 0;

  int x;
  __m128 a, b, c, d, e, f;
  __m128 three = _mm_set_ps1(3);

  b = _mm_load_ps((float*)&in[read++]);
  a = _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 0, 0, 0));
  c = _mm_load_ps((float*)&in[read++]);
  d = _mm_load_ps((float*)&in[read++]);

  for (x = 0; x < m128_pitch_out; ++x) {
    e = _mm_shuffle_ps(a, c, _MM_SHUFFLE(0, 0, 3, 3));
    f = _mm_shuffle_ps(b, d, _MM_SHUFFLE(0, 0, 3, 3));
    e = _mm_blend_ps(b, e, 9);
    f = _mm_blend_ps(c, f, 9);
    e = _mm_shuffle_ps(e, e, _MM_SHUFFLE(3, 1, 2, 0));
    f = _mm_shuffle_ps(f, f, _MM_SHUFFLE(3, 1, 2, 0));
    e = _mm_hadd_ps(e, f);

    f = _mm_hadd_ps(b, c);
    if (mul) {
      _mm_store_ps((float*)&Out[x],
                   _mm_mul_ps(_mm_add_ps(e, _mm_mul_ps(f, three)), *mul));
    } else {
      _mm_store_ps((float*)&Out[x], _mm_add_ps(e, _mm_mul_ps(f, three)));
    }
    a = c;
    b = d;
    if (read < m128_pitch_in - 1) {
      c = _mm_load_ps((float*)&in[read++]);
      d = _mm_load_ps((float*)&in[read++]);
    } else {
      if (read < m128_pitch_in) {
        c = _mm_load_ps((float*)&in[read++]);
        d = _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 3, 3, 3));
      } else {
        c = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 3, 3, 3));
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
  __m128* hi;
  __m128* lo;
  const __m128 four = _mm_set_ps1(4);
  const __m128 six = _mm_set_ps1(6);
  const __m128 eleven = _mm_set_ps1(11);
  const __m128 fifteen = _mm_set_ps1(15);
  const __m128 _16th = _mm_set_ps1(1.0 / 16);
  const __m128 _256th = _mm_set_ps1(1.0 / 256);

  for (l = 0; l < (int)levels_.size() - 1; ++l) {
    hi = (__m128*)levels_[l].data;
    lo = (__m128*)levels_[l + 1].data;

    int height_odd = (levels_[l].height & 1) ^ levels_[l].y_shift;
    int first_bad_line = levels_[l + 1].height - (3 - height_odd);

    for (int t = 0; t < (int)levels_[l + 1].bands.size() - 1; ++t) {
      threadpool_->Queue([=, this] {
        ShrinkThread(lines_[t], hi, lo, levels_[l].m128_pitch,
                     levels_[l + 1].m128_pitch, first_bad_line, height_odd,
                     levels_[l + 1].bands[t], levels_[l + 1].bands[t + 1],
                     levels_[l].x_shift, levels_[l].y_shift);
      });
    }

    threadpool_->Wait();
  }
}

void Pyramid::ShrinkThread(__m128* line, __m128* hi, __m128* lo,
                           int m128_pitch_hi, int m128_pitch_lo,
                           int first_bad_line, int height_odd, int sy, int ey,
                           const bool x_shift, const bool y_shift) {
  int x, y;

  const __m128 four = _mm_set_ps1(4);
  const __m128 six = _mm_set_ps1(6);
  const __m128 eleven = _mm_set_ps1(11);
  const __m128 fifteen = _mm_set_ps1(15);
  const __m128 _16th = _mm_set_ps1(1.0 / 16);
  const __m128 _256th = _mm_set_ps1(1.0 / 256);

  // line 0
  if (sy == 0) {
    memcpy(line, hi, m128_pitch_hi << 4);
    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _16th, x_shift);
    lo += m128_pitch_lo;

    if (!y_shift) {
      // line 1
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] =
            _mm_add_ps(_mm_add_ps(_mm_mul_ps(hi[x], eleven),
                                  _mm_mul_ps(hi[x + m128_pitch_hi], four)),
                       hi[x + (m128_pitch_hi << 1)]);
      }
      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      lo += m128_pitch_lo;
      hi += (m128_pitch_hi << 1);
      sy = 2;
    } else {
      // line 1
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] = _mm_add_ps(_mm_mul_ps(hi[x], fifteen), hi[x + m128_pitch_hi]);
      }
      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      lo += m128_pitch_lo;
      hi += m128_pitch_hi;
      // line 2
      for (x = 0; x < m128_pitch_hi; ++x) {
        _mm_store_ps(
            (float*)&line[x],
            _mm_add_ps(_mm_add_ps(_mm_add_ps(hi[x - m128_pitch_hi],
                                             hi[x + (m128_pitch_hi << 1)]),
                                  _mm_mul_ps(_mm_add_ps(hi[x - m128_pitch_hi],
                                                        hi[x + m128_pitch_hi]),
                                             four)),
                       _mm_mul_ps(hi[x], six)));
      }
      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      lo += m128_pitch_lo;
      hi += (m128_pitch_hi << 1);
      sy = 3;
    }
  } else {
    lo += m128_pitch_lo * sy;
    hi += m128_pitch_hi *
          (((sy - 1) << 1) -
           y_shift);  // because hi has a missing line compared to lo, lo's line
                      // 0 corresponds to hi's line -2, 1=>0, 2=>2, 3=>4
  }

  ey = (std::min)(first_bad_line, ey);

  // good lines
  for (y = sy; y < ey; ++y) {  // was y < first_bad_line
    for (x = 0; x < m128_pitch_hi; ++x) {
      _mm_store_ps(
          (float*)&line[x],
          _mm_add_ps(_mm_add_ps(_mm_add_ps(hi[x - (m128_pitch_hi << 1)],
                                           hi[x + (m128_pitch_hi << 1)]),
                                _mm_mul_ps(_mm_add_ps(hi[x - m128_pitch_hi],
                                                      hi[x + m128_pitch_hi]),
                                           four)),
                     _mm_mul_ps(hi[x], six)));
    }
    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
    lo += m128_pitch_lo;
    hi += (m128_pitch_hi << 1);
  }

  if (y == first_bad_line) {  // final block
    // prepenultimate line
    if (!height_odd) {
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] = _mm_add_ps(
            _mm_add_ps(
                _mm_add_ps(hi[x - (m128_pitch_hi << 1)], hi[x + m128_pitch_hi]),
                _mm_mul_ps(
                    _mm_add_ps(hi[x - m128_pitch_hi], hi[x + m128_pitch_hi]),
                    four)),
            _mm_mul_ps(hi[x], six));
      }

      Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
      ++y;
      lo += m128_pitch_lo;
      hi += (m128_pitch_hi << 1);

      // this case moved from block below
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] = _mm_add_ps(hi[x - (m128_pitch_hi << 1)],
                             _mm_mul_ps(hi[x - m128_pitch_hi], fifteen));
      }
    } else {
      // penultimate line
      for (x = 0; x < m128_pitch_hi; ++x) {
        line[x] =
            _mm_add_ps(_mm_add_ps(_mm_mul_ps(hi[x], eleven),
                                  _mm_mul_ps(hi[x - m128_pitch_hi], four)),
                       hi[x - (m128_pitch_hi << 1)]);
      }
      // other case removed from here, moved into block above
    }

    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _256th, x_shift);
    ++y;
    lo += m128_pitch_lo;
    hi += (m128_pitch_hi << 1);

    // last line
    hi -= m128_pitch_hi * (3 - height_odd);
    memcpy(line, hi, m128_pitch_hi << 4);
    Squeeze(line, lo, m128_pitch_lo, m128_pitch_hi, _16th, x_shift);
    ++y;
  }
}

void Pyramid::Squeeze(__m128* line, __m128* lo, int m128_pitch_lo,
                      int m128_pitch_hi, __m128 final_mul, bool x_shift) {
  int hi_x = 0;
  int lo_x = 0;

  __m128 a, b, c, d, e, f, g, h, i, j;

  const __m128 four = _mm_set_ps1(4);
  const __m128 six = _mm_set_ps1(6);

  if (x_shift) memmove(&((float*)line)[1], line, (m128_pitch_hi << 4) - 4);

  while (lo_x < m128_pitch_lo) {
    if (hi_x >= m128_pitch_hi) {  // was >= ... + 1
      b = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 3, 3, 3));
      c = b;
    } else {
      b = line[hi_x++];
      c = line[hi_x++];
      if (lo_x == 0) a = _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 0, 0, 0));
    }

    // a = EFGH
    // b = IJKL
    // c = MNOP

    // shuffle to four pairs of outer pixels
    f = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 0, 2, 0));  // EGIK
    g = _mm_shuffle_ps(b, c, _MM_SHUFFLE(2, 0, 2, 0));  // IKMO
    d = _mm_shuffle_ps(f, f, _MM_SHUFFLE(3, 1, 2, 0));  // EIGK
    e = _mm_shuffle_ps(g, g, _MM_SHUFFLE(3, 1, 2, 0));  // IMKO
    d = _mm_hadd_ps(d, e);

    // shuffle to four pairs of inner pixels
    h = _mm_shuffle_ps(
        a, b, _MM_SHUFFLE(1, 1, 3, 1));  // FHJJ // prev BDFF (should be BDDF)
    i = _mm_shuffle_ps(
        b, c, _MM_SHUFFLE(1, 1, 3, 1));  // JLNN // prev FHJJ (should be FHHJ)
    h = _mm_shuffle_ps(h, h, _MM_SHUFFLE(3, 1, 1, 0));  // FHHJ
    i = _mm_shuffle_ps(i, i, _MM_SHUFFLE(3, 1, 1, 0));  // JLLN
    h = _mm_mul_ps(four, _mm_hadd_ps(h, i));

    // shuffle to four central pixels
    j = _mm_mul_ps(six, _mm_shuffle_ps(f, g, _MM_SHUFFLE(2, 1, 2, 1)));  // GIKM

    // store
    _mm_store_ps((float*)&lo[lo_x++],
                 _mm_mul_ps(_mm_add_ps(d, _mm_add_ps(h, j)), final_mul));
    a = c;
  }
}

/***********************************************************************
 * laplace (gaussian)
 ***********************************************************************/
void Pyramid::LaplaceCollapse(int n_levels, bool Collapse) {
  int j, l;

  for (j = 0; j < n_levels - 1; ++j) {
    if (Collapse)
      l = (n_levels - 2) - j;
    else
      l = j;

    for (int t = 0; t < (int)levels_[l].bands.size() - 1; ++t) {
      threadpool_->Queue([=, this] {
        LaplaceThreadWrapper(&levels_[l], &levels_[l + 1], levels_[l].bands[t],
                             levels_[l].bands[t + 1]);
      });
    }

    threadpool_->Wait();
  }
}

void Pyramid::LaplaceExpand(__m128* hi, __m128* lo, int m128_pitch_hi,
                            int m128_pitch_lo) {
  __m128 p, q;

  const __m128 expand0 = _mm_set_ps(0, 0.125f, 0.75f, 0.125f);
  const __m128 expand1 = _mm_set_ps(0, 0.5f, 0.5f, 0);
  const __m128 expand2 = _mm_set_ps(0.125, 0.75f, 0.125f, 0);
  const __m128 expand3 = _mm_set_ps(0.5f, 0.5f, 0, 0);

  int x_hi = 0;
  int x_lo = 0;

  p = _mm_load_ps((float*)lo);
  ++x_lo;

  while (x_hi < m128_pitch_hi) {
    _mm_store_ps(
        (float*)&hi[x_hi],
        _mm_hadd_ps(
            _mm_hadd_ps(_mm_mul_ps(p, expand0), _mm_mul_ps(p, expand1)),
            _mm_hadd_ps(_mm_mul_ps(p, expand2), _mm_mul_ps(p, expand3))));

    ++x_hi;

    if (x_lo < m128_pitch_lo) {
      q = lo[x_lo];  // _mm_load_ps((float*)&lo[x_lo]);
    } else if (x_lo == m128_pitch_lo) {
      q = _mm_shuffle_ps(q, q, _MM_SHUFFLE(3, 3, 3, 3));
    }
    ++x_lo;
    p = _mm_shuffle_ps(p, q, _MM_SHUFFLE(1, 0, 3, 2));

    _mm_store_ps(
        (float*)&hi[x_hi],
        _mm_hadd_ps(
            _mm_hadd_ps(_mm_mul_ps(p, expand0), _mm_mul_ps(p, expand1)),
            _mm_hadd_ps(_mm_mul_ps(p, expand2), _mm_mul_ps(p, expand3))));

    p = q;

    ++x_hi;
  }
}

void Pyramid::LaplaceExpandShifted(__m128* hi, __m128* lo, int m128_pitch_hi,
                                   int m128_pitch_lo) {
  __m128 p, q, t;

  const __m128 expand0 = _mm_set_ps(0, 0.5f, 0.5f, 0);
  const __m128 expand1 = _mm_set_ps(0.125, 0.75f, 0.125f, 0);
  const __m128 expand2 = _mm_set_ps(0.5f, 0.5f, 0, 0);
  const __m128 expand3 = _mm_set_ps(0, 0.125f, 0.75f, 0.125f);

  int x_hi = 0;
  int x_lo = 0;

  p = _mm_load_ps((float*)lo);
  ++x_lo;

  while (x_hi < m128_pitch_hi) {
    t = p;

    if (x_lo < m128_pitch_lo) {
      q = lo[x_lo];
    } else if (x_lo == m128_pitch_lo) {
      q = _mm_shuffle_ps(q, q, _MM_SHUFFLE(3, 3, 3, 3));
    }
    ++x_lo;
    p = _mm_shuffle_ps(p, q, _MM_SHUFFLE(1, 0, 3, 2));

    _mm_store_ps(
        (float*)&hi[x_hi],
        _mm_hadd_ps(
            _mm_hadd_ps(_mm_mul_ps(t, expand0), _mm_mul_ps(t, expand1)),
            _mm_hadd_ps(_mm_mul_ps(t, expand2), _mm_mul_ps(p, expand3))));

    ++x_hi;

    _mm_store_ps(
        (float*)&hi[x_hi],
        _mm_hadd_ps(
            _mm_hadd_ps(_mm_mul_ps(p, expand0), _mm_mul_ps(p, expand1)),
            _mm_hadd_ps(_mm_mul_ps(p, expand2), _mm_mul_ps(q, expand3))));

    p = q;

    ++x_hi;
  }
}

void Pyramid::LaplaceLine2(__m128* hi, __m128* temp1, __m128* temp2,
                           int m128_pitch) {
  const __m128 half = _mm_set_ps1(0.5f);

  for (int x = 0; x < m128_pitch; ++x) {
    _mm_store_ps(
        (float*)&hi[x],
        _mm_sub_ps(_mm_mul_ps(_mm_add_ps(temp1[x], temp2[x]), half), hi[x]));
  }
}

void Pyramid::LaplaceLine3(__m128* hi, __m128* temp1, __m128* temp2,
                           __m128* temp3, int m128_pitch) {
  static const __m128 eighth = _mm_set_ps1(0.125f);
  static const __m128 three_quarters = _mm_set_ps1(0.75f);

  for (int x = 0; x < m128_pitch; ++x) {
    _mm_store_ps(
        (float*)&hi[x],
        _mm_sub_ps(
            _mm_add_ps(_mm_mul_ps(_mm_add_ps(temp1[x], temp3[x]), eighth),
                       _mm_mul_ps(temp2[x], three_quarters)),
            hi[x]));
  }
}

__m128* GetLine(const Pyramid::Level& level, int y) {
  return (__m128*)(level.data + y * level.pitch);
}

void GetExpandedLine(const Pyramid::Level& level, __m128* temp, int y) {
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

  __m128* temp1 = (__m128*)_aligned_malloc(temp, 16);
  __m128* temp2 = (__m128*)_aligned_malloc(temp, 16);
  __m128* temp3 = (__m128*)_aligned_malloc(temp, 16);

  LaplaceThread(upper_level, lower_level, sy, ey, temp1, temp2, temp3);

  _aligned_free(temp1);
  _aligned_free(temp2);
  _aligned_free(temp3);
}

void Pyramid::LaplaceThread(Level* upper_level, Level* lower_level, int sy,
                            int ey, __m128* temp1, __m128* temp2,
                            __m128* temp3) {
  __m128* hi = (__m128*)upper_level->data + sy * upper_level->m128_pitch;

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
    if (!((y + upper_level->y_shift) & 1)) {
      GetExpandedLine(*lower_level, temp3, lo_y++);

      LaplaceLine3(hi, temp1, temp2, temp3, upper_level->m128_pitch);
    } else {
      LaplaceLine2(hi, temp2, temp3, upper_level->m128_pitch);

      __m128* temp = temp1;
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
  int x, y;
  int fours = levels_[0].width >> 2;

  __m128 m128_total = {0};
  __m128 one = _mm_set_ps1(1.0f);
  double total = 0;
  double row_total;

  __m128* data = (__m128*)levels_[0].data;

  for (y = 0; y < levels_[0].height; ++y) {
    m128_total = _mm_setzero_ps();

    for (x = 0; x < fours; ++x) {
      m128_total = _mm_add_ps(m128_total, data[x]);
    }

    m128_total = _mm_hadd_ps(m128_total, m128_total);
    m128_total = _mm_hadd_ps(m128_total, m128_total);
    row_total = _mm_cvtss_f32(m128_total);

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
  __m128 __add = _mm_set_ps1(add);

  int lim = (std::min)(levels, (int)levels_.size() - 1);

  for (int l = 0; l < lim; ++l) {
    __m128* data = (__m128*)levels_[l].data;

    for (int t = 0; t < (int)levels_[l].bands.size() - 1; ++t) {
      threadpool_->Queue([=, this]() {
        __m128* data = (__m128*)levels_[l].data +
                       levels_[l].bands[t] * levels_[l].m128_pitch;
        for (int y = levels_[l].bands[t]; y < levels_[l].bands[t + 1]; ++y) {
          for (int x = 0; x < levels_[l].m128_pitch; ++x) {
            _mm_store_ps((float*)&data[x], _mm_add_ps(__add, data[x]));
          }
          data += levels_[l].m128_pitch;
        }
      });
    }

    threadpool_->Wait();
  }
}

/***********************************************************************
 * Multiply and add
 ***********************************************************************/
void Pyramid::MultiplyAndAdd(float add, float mul, int levels) {
  int i;
  int x, y;
  __m128 __add = _mm_set_ps1(add);
  __m128 __mul = _mm_set_ps1(mul);

  int lim = (std::min)(levels, (int)levels_.size() - 1);

  for (i = 0; i < lim; ++i) {
    __m128* data = (__m128*)levels_[i].data;

    for (y = 0; y < levels_[i].height; ++y) {
      for (x = 0; x < levels_[i].m128_pitch; ++x) {
        _mm_store_ps((float*)&data[x],
                     _mm_add_ps(__add, _mm_mul_ps(data[x], __mul)));
      }
      data += levels_[i].m128_pitch;
    }
  }
}

/***********************************************************************
 * Multiply, add, clamp
 ***********************************************************************/
void Pyramid::MultiplyAddClamp(float add, float mul, int level) {
  int x, y;
  __m128 __add = _mm_set_ps1(add);
  __m128 __mul = _mm_set_ps1(mul);
  __m128 __min = _mm_set_ps1(1.0f);
  __m128 __max = _mm_set_ps1(0.0f);

  __m128* data = (__m128*)levels_[level].data;

  for (y = 0; y < levels_[level].height; ++y) {
    for (x = 0; x < levels_[level].m128_pitch; ++x) {
      _mm_store_ps(
          (float*)&data[x],
          _mm_max_ps(
              _mm_min_ps(_mm_add_ps(__add, _mm_mul_ps(data[x], __mul)), __min),
              __max));
    }
    data += levels_[level].m128_pitch;
  }
}

/***********************************************************************
 * Multiply
 ***********************************************************************/
void Pyramid::Multiply(int level, float mul) {
  if (mul == 1) return;
  if (mul == 0) {
    ZeroMemory(levels_[level].data,
               levels_[level].height * levels_[level].pitch * sizeof(float));
    return;
  }

  int x, y;
  __m128 __mul = _mm_set_ps1(mul);

  __m128* data = (__m128*)levels_[level].data;

  for (y = 0; y < levels_[level].height; ++y) {
    for (x = 0; x < levels_[level].m128_pitch; ++x) {
      _mm_store_ps((float*)&data[x],
                   _mm_mul_ps(_mm_load_ps((float*)&data[x]), __mul));
    }
    data += levels_[level].m128_pitch;
  }
}

/***********************************************************************
 * multiply_by_pyramid
 ***********************************************************************/
void Pyramid::MultplyByPyramid(Pyramid* b) {
  int x, y;

  for (int l = 0; l < (int)levels_.size() - 1; ++l) {
    __m128* data = (__m128*)levels_[l].data;
    __m128* _b = (__m128*)b->levels_[l].data;

    for (y = 0; y < levels_[l].height; ++y) {
      for (x = 0; x < levels_[l].m128_pitch; ++x) {
        _mm_store_ps((float*)&data[x], _mm_mul_ps(_mm_load_ps((float*)&data[x]),
                                                  _mm_load_ps((float*)&_b[x])));
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
    //  fuse_thread((__m128*)levels[l].data,
    //(__m128*)_b->levels[l].data,
    //(__m128*)mask->levels[l].data, m128_pitch, 0, levels[l].height, pre,
    // black);

    // fuse doesn't see any gains from multithreading; leave this here as
    // reference

    for (int t = 0; t < (int)levels_[l].bands.size() - 1; ++t) {
      threadpool_->Queue([=, this] {
        FuseThread((__m128*)levels_[l].data, (__m128*)_b->levels_[l].data,
                   (__m128*)mask->levels_[l].data, levels_[l].m128_pitch,
                   levels_[l].bands[t], levels_[l].bands[t + 1], pre, black);
      });
    }
    threadpool_->Wait();
  }
}

void Pyramid::FuseThread(__m128* a, __m128* b, __m128* m, int m128_pitch,
                         int sy, int ey, bool pre, int black) {
  int p;
  int add = sy * m128_pitch;
  int count = (ey - sy) * m128_pitch;

  a += add;
  b += add;
  m += add;

  if (!pre) {
    for (p = 0; p < count; ++p) {
      __m128 _a = a[p];
      _mm_store_ps((float*)&a[p],
                   _mm_add_ps(_a, _mm_mul_ps(_mm_sub_ps(b[p], _a), m[p])));
    }
  } else {
    __m128 ones = _mm_set_ps1(1.0f);
    __m128 blacks = _mm_set_ps1((float)black);
    if (black) {
      for (p = 0; p < count; ++p) {
        _mm_store_ps(
            (float*)&a[p],
            _mm_add_ps(blacks, _mm_add_ps(_mm_sub_ps(b[p], blacks),
                                          _mm_mul_ps(_mm_sub_ps(a[p], blacks),
                                                     _mm_sub_ps(ones, m[p])))));
      }
    } else {
      for (p = 0; p < count; ++p) {
        _mm_store_ps(
            (float*)&a[p],
            _mm_add_ps(b[p], _mm_mul_ps(a[p], _mm_sub_ps(ones, m[p]))));
      }
    }
  }
}

void Pyramid::Fuse(Pyramid* b, float weight) {
  int l;
  int p;
  __m128 w = _mm_set_ps1(weight);

  for (l = 0; l < (int)levels_.size(); ++l) {
    __m128* _a = (__m128*)levels_[l].data;
    __m128* _b = (__m128*)b->levels_[l].data;

    int count = levels_[l].height * levels_[l].m128_pitch;
    for (p = 0; p < count; ++p) {
      __m128 __a = _a[p];
      _mm_store_ps((float*)&_a[p],
                   _mm_add_ps(__a, _mm_mul_ps(_mm_sub_ps(_b[p], __a), w)));
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
  __m128 one = _mm_set_ps1(1);
  __m128 half = _mm_set_ps1(0.5f);
  __m128 _power = _mm_set_ps1(power);
  __m128 pi = _mm_set_ps1(3.14159265359f);
  __m128i andi = _mm_set_epi32(0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff);
  __m128* _and = (__m128*)&andi;

  if (gamma) _power = _mm_mul_ps(_power, _power);
  _power = _mm_div_ps(one, _power);

  __m128* data = (__m128*)levels[level].data;
  for (y = 0; y < levels[level].height; ++y) {
    for (x = 0; x < levels[level].m128_pitch; ++x) {
      __m128 d = data[x];
      _mm_store_ps(
          (float*)&data[x],
          _mm_mul_ps(
              _mm_mul_ps(
                  _mm_sub_ps(
                      one, cos_ps(_mm_min_ps(
                               _mm_and_ps(_mm_mul_ps(d, _power), *_and), pi))),
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
  if (b->GetNLevels() < GetNLevels()) return;
  memcpy(levels_[GetNLevels() - 1].data, b->levels_[GetNLevels() - 1].data,
         levels_[GetNLevels() - 1].height * levels_[GetNLevels() - 1].pitch *
             sizeof(float));
}

/***********************************************************************
 * approximate gaussian blur
 ***********************************************************************/
#define BLUR_SSE_GET(y, x) \
  y = _mm_set_ps(line3[x], line2[x], line1[x], line0[x])
// #define BLUR_SSE_GET2(x) _mm_load_ps((float*)&transposed[x])
#define BLUR_SSE_GET_LEFT                                                 \
  temp1 = _mm_set_ps(line3[left], line2[left], line1[left], line0[left]); \
  left++;
#define BLUR_SSE_GET_RIGHT                                                    \
  temp2 = _mm_set_ps(line3[right], line2[right], line1[right], line0[right]); \
  right++;

void Pyramid::BlurX(float radius, Pyramid* transpose) {
  for (int i = 0; i < (int)levels_[0].bands.size() - 1; ++i) {
    threadpool_->Queue([=, this] {
      BlurXThread(radius, transpose, levels_[0].bands[i],
                  levels_[0].bands[i + 1]);
    });
  }
  threadpool_->Wait();
}

void Pyramid::BlurXThread(float radius, Pyramid* transpose, int sy, int ey) {
  int x, y;
  int i;
  int o;
  float* line0 = levels_[0].data + sy * levels_[0].pitch;
  float* line1 = line0 + levels_[0].pitch;
  float* line2 = line1 + levels_[0].pitch;
  float* line3 = line2 + levels_[0].pitch;
  float* out = transpose->levels_[0].data + sy;
  __m128 temp1, temp2;

  int iradius = (int)floor(radius);
  __m128 irp1 = _mm_set_ps1((float)(iradius + 1));
  __m128 mul = _mm_set_ps1(radius - iradius);
  __m128 acc;

  int left, right;

  int fours = (ey - sy + 3) >>
              2;  // +3 is probably not necessary because all bands are mod 4

  if (iradius < levels_[0].width >> 1) {
    for (y = 0; y < fours; ++y) {
      acc = _mm_setzero_ps();
      left = 0;

      BLUR_SSE_GET_LEFT;

      acc = _mm_mul_ps(temp1, irp1);
      for (right = 1; right < iradius + 1;) {
        BLUR_SSE_GET_RIGHT;
        acc = _mm_add_ps(acc, temp2);
      }

      x = 0;
      o = 0;
      right = iradius + 1;

      for (i = 0; i <= iradius; ++i) {
        BLUR_SSE_GET_RIGHT;
        _mm_store_ps(
            &out[o],
            _mm_add_ps(acc, _mm_mul_ps(_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        acc = _mm_add_ps(_mm_sub_ps(temp2, temp1), acc);
        ++x;
      }

      while (right < levels_[0].width) {
        BLUR_SSE_GET_RIGHT;
        _mm_store_ps(
            &out[o],
            _mm_add_ps(acc, _mm_mul_ps(_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        BLUR_SSE_GET_LEFT
        acc = _mm_add_ps(_mm_sub_ps(temp2, temp1), acc);
        ++x;
      }

      while (x < levels_[0].width) {
        _mm_store_ps(
            &out[o],
            _mm_add_ps(acc, _mm_mul_ps(_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        BLUR_SSE_GET_LEFT;
        acc = _mm_add_ps(_mm_sub_ps(temp2, temp1), acc);
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
      acc = _mm_setzero_ps();

      BLUR_SSE_GET(temp1, 0);
      acc = _mm_mul_ps(temp1, irp1);
      right = 1;
      for (x = 1; x < iradius + 1; ++x) {
        if (right < levels_[0].width) {
          BLUR_SSE_GET(temp2, right);
          ++right;
        }
        acc = _mm_add_ps(acc, temp2);
      }

      x = 0;
      o = 0;
      left = -iradius;

      for (x = 0; x < levels_[0].width; ++x) {
        if (right < levels_[0].width) {
          BLUR_SSE_GET(temp2, right);
          ++right;
        }
        _mm_store_ps(
            &out[o],
            _mm_add_ps(acc, _mm_mul_ps(_mm_add_ps(temp1, temp2), mul)));
        o += transpose->levels_[0].pitch;
        if (left > 0) BLUR_SSE_GET(temp1, left);
        ++left;
        acc = _mm_add_ps(_mm_sub_ps(temp2, temp1), acc);
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
  uint8_t* temp = (uint8_t*)calloc(width * height, 1);

  int px = 0, py = 0;

  for (int l = 0; l < (int)levels_.size(); ++l) {
    float* data = (float*)levels_[l].data;
    uint8_t* line = temp + py * levels_[0].pitch + px;
    for (int y = 0; y < levels_[l].height; ++y) {
      for (int x = 0; x < levels_[l].pitch; ++x) {
        int f = (int)floor(data[x] + 0.5);
        line[x] = std::max(0, std::min(255, f));
      }
      line += levels_[0].pitch;
      data += levels_[l].pitch;
    }
    if (l & 1)
      px += levels_[l].pitch + 1;
    else
      py += levels_[l].height + 1;
  }

  io::png::Pnger::Quick((char*)filename, temp, width, height, width,
                        PNG_COLOR_TYPE_GRAY);

  free(temp);
}

}  // namespace multiblend
