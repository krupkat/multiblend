// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mb/functions.h"

#include <chrono>
#include <cstdarg>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

#include "mb/flex.h"
#include "mb/logging.h"
#include "mb/pyramid.h"

namespace multiblend::utils {

/***********************************************************************
 * ShrinkMasks
 ***********************************************************************/
int Squish(const uint32_t* in, uint32_t* out, int in_width, int out_width) {
  float current_val;
  uint32_t cur;
  float a;
  float b;
  float c;
  float d;
  float e;
  int last_int = -1;

  int in_p = 0;
  int out_p = 0;
  int read = 0;
  int wrote = 0;
  int in_count;
  int out_count;

  cur = in[in_p++];
  if ((cur & 0x80000000) != 0u) {
    in_count = cur & 0x00ffffff;
    if ((cur & 0x20000000) != 0u) {
      current_val = ((float*)in)[in_p++];
    } else {
      last_int = (cur >> 30) & 1;
    }
  } else {
    in_count = 1;
    current_val = *((float*)&cur);
  }

  d = -1;
  e = -1;

  // trivial: full width block
  if (in_count == in_width) {
    if (last_int >= 0) {
      *out = (cur & 0xff000000) | out_width;
    } else {
      out[out_p++] = (cur & 0xa0000000) | out_width;
      ((float*)out)[out_p] = current_val;
    }
    return in_p;
  }

  out_count = (in_count + 1) >> 1;
  if (last_int >= 0) {
    out[out_p++] = (cur & 0xff000000) | out_count;
    current_val = (float)last_int;
  } else {
    if (out_count > 1) {
      out[out_p++] = (cur & 0xa0000000) | out_count;
    }
    ((float*)out)[out_p++] = current_val;
  }
  wrote = out_count;

  a = b = c = current_val;
  read = in_count;
  in_count -= (in_count - 1) | 1;

  int c_read = in_p;
  int e_read;

  while (wrote < out_width) {
    // loop always starts needing to read two more pixels

    // read d
    if ((in_count == 0) && read < in_width) {
      cur = in[in_p++];
      if ((cur & 0x80000000) != 0u) {  // float
        in_count = cur & 0x00ffffff;
        read += in_count;
        if ((cur & 0x20000000) != 0u) {  // repeated float
          last_int = -1;
          current_val = *((float*)&in[in_p++]);
        } else {
          last_int = (cur >> 30) & 1;
          current_val = (float)last_int;
        }
      } else {
        current_val = *((float*)&cur);
        last_int = -1;
        in_count = 1;
        ++read;
      }
    }

    if (read == in_width) {
      in_count = 0x7fffffff;
    }

    if (in_count >= 2) {
      d = e = current_val;
      e_read = in_p;
      in_count -= 2;
    } else {
      d = current_val;

      cur = in[in_p++];
      if ((cur & 0x80000000) != 0u) {
        in_count = cur & 0x00ffffff;
        read += in_count;
        if ((cur & 0x20000000) != 0u) {
          last_int = -1;
          current_val = ((float*)in)[in_p++];
        } else {
          last_int = (cur >> 30) & 1;
          current_val = (float)last_int;
        }
      } else {
        current_val = *((float*)&cur);
        last_int = -1;
        in_count = 1;
        ++read;
      }

      e = current_val;
      e_read = in_p;
      --in_count;
    }

    // calculate
    ((float*)out)[out_p++] = (b + d) * 0.25f + (a + e + 6 * c) * 0.0625f;
    wrote++;

    if (c_read == in_p && in_count >= 4) {
      out_count = in_count >> 1;
      if (out_count > (out_width - wrote)) {
        out_count = out_width - wrote;
      }
      if (last_int >= 0) {
        out[out_p++] = ((0x2 | last_int) << 30) | out_count;
      } else {
        out[out_p++] = 0xa0000000 | out_count;
        ((float*)out)[out_p++] = e;
      }
      wrote += out_count;
      in_count -= out_count << 1;
    }
    a = c;
    b = d;
    c = e;
    c_read = e_read;
  }

  return in_p;
}

void ShrinkMasks(std::vector<Flex>& masks, int n_levels) {
  int i;
  Flex flex_temp(masks[0].width_, masks[0].height_);
  uint32_t cur;
  uint32_t* lines[5];
  std::unique_ptr<uint32_t[]> real_lines[5];
  int count[5];
  int pointer[5];
  float vals[5];
  int min_count;

  for (auto& real_line : real_lines) {
    real_line = std::make_unique<uint32_t[]>(masks[0].width_);
  }

  for (int l = 1; l < n_levels; ++l) {
    int in_width = masks[l - 1].width_;
    int out_width = (in_width + 6) >> 1;
    int out_height = (masks[l - 1].height_ + 6) >> 1;
    masks.emplace_back(out_width, out_height);

    int input_p = 0;

    // first line
    input_p += Squish((uint32_t*)&masks[l - 1].data_.get()[input_p],
                      real_lines[0].get(), in_width, out_width)
               << 2;
    int lines_read = 1;

    int wrote = 0;
    pointer[0] = 0;
    while (wrote < out_width) {
      cur = real_lines[0][pointer[0]++];
      masks[l].Write32(cur);
      if ((cur & 0x80000000) == 0u) {
        ++wrote;
      } else {
        if ((cur & 0x20000000) != 0u) {
          masks[l].Write32(real_lines[0][pointer[0]++]);
        }
        wrote += cur & 0x00ffffff;
      }
    }

    // other lines
    lines[0] = lines[1] = lines[2] = real_lines[0].get();
    lines[3] = real_lines[3].get();
    lines[4] = real_lines[4].get();

    input_p += Squish((uint32_t*)&masks[l - 1].data_.get()[input_p], lines[3],
                      in_width, out_width)
               << 2;
    input_p += Squish((uint32_t*)&masks[l - 1].data_.get()[input_p], lines[4],
                      in_width, out_width)
               << 2;
    lines_read += 2;

    for (int y = 1; y < masks[l].height_; ++y) {
      masks[l].NextLine();

      wrote = 0;
      pointer[0] = pointer[1] = pointer[2] = pointer[3] = pointer[4] = 0;
      count[0] = count[1] = count[2] = count[3] = count[4] = 0;
      int last_int;

      while (wrote < out_width) {
        min_count = 0x7fffffff;
        for (i = 0; i < 5; ++i) {
          if (count[i] == 0) {
            cur = lines[i][pointer[i]++];
            if ((cur & 0x80000000) == 0u) {
              vals[i] = *((float*)&cur);
              count[i] = 1;
            } else {
              if ((cur & 0x20000000) != 0u) {
                vals[i] = ((float*)lines[i])[pointer[i]++];
              } else {
                last_int = (cur >> 30) & 1;
                vals[i] = (float)last_int;
              }
              count[i] = cur & 0x00ffffff;
            }
          }

          if (count[i] < min_count) {
            min_count = count[i];
          }
        }

        float val = (vals[0] + vals[4]) * 0.0625f +
                    (vals[1] + vals[3]) * 0.25f + vals[2] * 0.375f;
        if (val == 0 || val == 1) {
          if (min_count == 1) {
            masks[l].Write32(*((uint32_t*)&val));
          } else {
            int _int = (int)val;
            masks[l].Write32(0x80000000 | (_int << 30) | min_count);
          }
        } else {
          if (min_count > 1) {
            masks[l].Write32(0xa0000000 | min_count);
          }
          masks[l].Write32(*((uint32_t*)&val));
        }

        wrote += min_count;
        for (i = 0; i < 5; ++i) {
          count[i] -= min_count;
        }
      }

      if (y == 1) {
        lines[0] = real_lines[1].get();
        lines[1] = real_lines[2].get();
      }

      uint32_t* temp = lines[0];
      lines[0] = lines[2];
      lines[2] = lines[4];
      lines[4] = lines[1];
      lines[1] = lines[3];
      lines[3] = temp;

      if (lines_read < masks[l - 1].height_) {
        input_p += Squish((uint32_t*)&masks[l - 1].data_.get()[input_p],
                          lines[3], in_width, out_width)
                   << 2;
        ++lines_read;
      } else {
        lines[3] = lines[2];
      }

      if (lines_read < masks[l - 1].height_) {
        input_p += Squish((uint32_t*)&masks[l - 1].data_.get()[input_p],
                          lines[4], in_width, out_width)
                   << 2;
        ++lines_read;
      } else {
        lines[4] = lines[3];
      }
    }
  }
}

/***********************************************************************
 * Composite line
 ***********************************************************************/
void CompositeLine(const float* input_p, float* output_p, int i, int x_offset,
                   int in_level_width, int out_level_width, int out_level_pitch,
                   uint8_t* _mask, std::size_t mask_p) {
  int x = 0;
  int mask_count;
  float current_val;
  int last_int;

  mask_p >>= 2;
  auto* mask = (uint32_t*)_mask;

  while (x < out_level_width) {
    uint32_t cur = mask[mask_p++];
    last_int = -1;
    if ((cur & 0x80000000) != 0u) {
      mask_count = cur & 0x0fffffff;
      if ((cur & 0x20000000) != 0u) {
        current_val = ((float*)mask)[mask_p++];
      } else {
        last_int = (cur >> 30) & 1;
      }
    } else {
      mask_count = 1;
      current_val = *((float*)&cur);
    }

    int lim = x + mask_count;

    if (last_int == 0) {
      if (i == 0) {
        memset(&output_p[x], 0, mask_count << 2);
      }
      x += mask_count;
    } else if (last_int == 1) {
      if (x < x_offset) {
        float f = input_p[0];
        while (x < lim && x < x_offset) {
          output_p[x++] = f;
        }
      }

      while (x < lim && x < x_offset + in_level_width) {
        output_p[x] = input_p[x - x_offset];
        ++x;
      }

      if (x < lim) {
        float f = input_p[in_level_width - 1];
        while (x < lim) {
          output_p[x++] = f;
        }
      }
    } else {
      if (x < x_offset) {
        float f = input_p[0] * current_val;
        while (x < lim && x < x_offset) {
          if (i == 0) {
            output_p[x++] = f;
          } else {
            output_p[x++] += f;
          }
        }
      }

      while (x < lim && x < x_offset + in_level_width) {
        float f = input_p[x - x_offset] * current_val;
        if (i == 0) {
          output_p[x++] = f;
        } else {
          output_p[x++] += f;
        }
      }

      if (x < lim) {
        float f = input_p[in_level_width - 1] * current_val;
        while (x < lim) {
          if (i == 0) {
            output_p[x++] = f;
          } else {
            output_p[x++] += f;
          }
        }
      }
    }
  }

  if (x < out_level_pitch) {
    float f = output_p[x - 1];
    do {
      output_p[x++] = f;
    } while (x < out_level_pitch);
  }
}

/***********************************************************************
 * DT reader and macros
 ***********************************************************************/
void ReadInpaintDT(Flex& flex, int& current_count, int& current_step,
                   uint32_t& dt_val) {
  if (current_count != 0) {
    --current_count;
    dt_val += current_step;
  } else {
    uint8_t _byte = flex.ReadBackwards8();
    if (_byte == 255) {
      dt_val = flex.ReadBackwards32();
      return;
    }
    current_step = ((_byte & 7) - 3);
    if ((_byte & 0x80) == 0) {  // 0b0ccccsss
      current_count = _byte >> 3;
    } else if ((_byte & 0x40) == 0) {  // 0b10ssssss
      current_step = (_byte & 0x3f);
      current_count = 0;
    } else if ((_byte & 0x20) == 0) {  // 0b11000000
      current_count = flex.ReadBackwards8();
    } else if ((_byte & 0x10) == 0) {  // 0b11100000
      current_count = flex.ReadBackwards16();
    } else {  // if (!(_byte & 0x08)) {
      current_count = flex.ReadBackwards32();
    }
    dt_val += current_step;
  }
}

void ReadSeamDT(Flex& flex, int& current_count, int64_t& current_step,
                uint64_t& dt_val) {
  if (current_count != 0) {
    --current_count;
    dt_val += current_step;
  } else {
    uint8_t _byte = flex.ReadBackwards8();
    if (_byte == 255) {
      dt_val = flex.ReadBackwards64();
      return;
    }
    current_step = ((int64_t)(_byte & 7) - 3) << 32;
    if ((_byte & 0x80) == 0) {  // 0b0ccccsss
      current_count = _byte >> 3;
    } else if ((_byte & 0x40) == 0) {  // 0b10ssssss
      current_step = (int64_t)(_byte & 0x3f) << 32;
      current_count = 0;
    } else if ((_byte & 0x20) == 0) {  // 0b11000000
      current_count = flex.ReadBackwards8();
    } else if ((_byte & 0x10) == 0) {  // 0b11100000
      current_count = flex.ReadBackwards16();
    } else {  // if (!(_byte & 0x08)) {
      current_count = flex.ReadBackwards32();
    }
    dt_val += current_step;
  }
}

/***********************************************************************
 * Seam line compress
 ***********************************************************************/
void cwrite(int current_count, int current_step, uint8_t* output, int& p) {
  if (current_count != 0) {
    if (--current_count < 16) {
      output[p++] = current_count << 3 | current_step;
    } else if (current_count < 256) {
      output[p++] = current_count;
      output[p++] = 0xc0 | current_step;
    } else if (current_count < 65536) {
      *((uint16_t*)&output[p]) = current_count;
      p += 2;
      output[p++] = 0xe0 | current_step;
    } else {
      *((uint32_t*)&output[p]) = current_count;
      p += 4;
      output[p++] = 0xf0 | current_step;
    }
  }
}

int CompressDTLine(const uint32_t* input, uint8_t* output, int width) {
  int current_step = -100;
  int current_count = 0;

  int step;
  int x = 0;
  int p = 0;
  uint32_t left_val;
  uint32_t right_val;

  while (((left_val = input[x++]) == 0u) && x < width) {
    ;
  }
  if (left_val == 0u) {
    return p;
  }

  while (x < width) {
    while (((right_val = input[x++]) == 0u) && x < width) {
      ;
    }
    if (right_val == 0u) {
      break;
    }

    if ((step = (int)(left_val - right_val) + 3) < 67 && step >= 0) {
      if (step <= 7) {
        if (step == current_step) {
          ++current_count;
        } else {
          cwrite(current_count, current_step, output, p);
          current_step = step;
          current_count = 1;
        }
      } else {
        cwrite(current_count, current_step, output, p);
        output[p++] = 0x80 | (step - 3);
        current_step = -100;
        current_count = 0;
      }
    } else {
      cwrite(current_count, current_step, output, p);
      *((uint32_t*)&output[p]) = left_val;
      p += 4;
      output[p++] = -1;
      current_step = -100;
      current_count = 0;
    }

    left_val = right_val;
  }

  cwrite(current_count, current_step, output, p);
  *((uint32_t*)&output[p]) = left_val;
  p += 4;

  output[p++] = -1;

  return p;
}

int CompressSeamLine(const uint64_t* input, uint8_t* output, int width) {
  int current_step = -100;
  int current_count = 0;

  int64_t step;
  int x = width;
  int p = 0;
  uint64_t left_val;
  uint64_t right_val;

  while ((((right_val = input[--x]) & 0xffffffff00000000) == 0u) && x > 0) {
    ;
  }
  if ((right_val & 0xffffffff00000000) == 0u) {
    return p;
  }

  while (x > 0) {
    while ((((left_val = input[--x]) & 0xffffffff00000000) == 0u) && x > 0) {
      ;
    }
    if ((left_val & 0xffffffff00000000) == 0u) {
      break;
    }

    if ((((right_val ^ left_val) & 0xffffffff) == 0u) &&
        (step = ((int64_t)(right_val - left_val) >> 32) + 3) < 67 &&
        step >= 0) {  // was <= 7
      if (step <= 7) {
        if (step == current_step) {
          ++current_count;
        } else {
          cwrite(current_count, current_step, output, p);
          current_step = (int)step;
          current_count = 1;
        }
      } else {
        cwrite(current_count, current_step, output, p);
        output[p++] = 0x80 | ((int)step - 3);
        current_step = -100;
        current_count = 0;
      }
    } else {
      cwrite(current_count, current_step, output, p);
      *((uint64_t*)&output[p]) = right_val;
      p += 8;
      output[p++] = -1;
      current_step = -100;
      current_count = 0;
    }

    right_val = left_val;
  }

  cwrite(current_count, current_step, output, p);
  *((uint64_t*)&output[p]) = right_val;
  p += 8;

  output[p++] = -1;

  return p;
}

/***********************************************************************
 * Wrap juggling
 ***********************************************************************/
void SwapUnswapH(Pyramid* py, bool unswap) {
  int width = py->GetWidth();
  int height = py->GetHeight();

  int minor_bytes = (width >> 1) << 2;
  int major_bytes = ((width + 1) >> 1) << 2;
  float* data = py->GetData();
  auto temp = std::make_unique<uint8_t[]>(major_bytes);

  for (int y = 0; y < height; ++y) {
    if (unswap) {
      memcpy(temp.get(), &((uint8_t*)data)[minor_bytes], major_bytes);
      memcpy(&((uint8_t*)data)[major_bytes], data, minor_bytes);
      memcpy(data, temp.get(), major_bytes);
    } else {
      memcpy(temp.get(), data, major_bytes);
      memcpy(data, &((uint8_t*)data)[major_bytes], minor_bytes);
      memcpy(&((uint8_t*)data)[minor_bytes], temp.get(), major_bytes);
    }

    data += py->GetPitch();
  }
}

void SwapH(Pyramid* py) { SwapUnswapH(py, false); }

void UnswapH(Pyramid* py) { SwapUnswapH(py, true); }

void SwapUnswapV(Pyramid* py, bool unswap) {
  int height = py->GetHeight();
  int byte_pitch = py->GetPitch() << 2;
  auto temp = std::make_unique<uint8_t[]>(byte_pitch);

  int half_height = height >> 1;

  if ((height & 1) != 0) {
    auto temp2 = std::make_unique<uint8_t[]>(byte_pitch);
    if (unswap) {
      auto* upper =
          (uint8_t*)(py->GetData() + ((int64_t)height >> 1) * py->GetPitch());
      auto* lower =
          (uint8_t*)(py->GetData() + ((int64_t)height - 1) * py->GetPitch());

      memcpy(temp.get(), upper, byte_pitch);

      for (int y = 0; y < half_height; ++y) {
        memcpy(temp2.get(), lower, byte_pitch);
        memcpy(lower, temp.get(), byte_pitch);
        upper -= byte_pitch;
        memcpy(temp.get(), upper, byte_pitch);
        memcpy(upper, temp2.get(), byte_pitch);
        lower -= byte_pitch;
      }

      memcpy(lower, temp.get(), byte_pitch);
    } else {
      auto* upper = (uint8_t*)py->GetData();
      auto* lower =
          (uint8_t*)(py->GetData() + ((int64_t)height >> 1) * py->GetPitch());

      memcpy(temp.get(), lower, byte_pitch);

      for (int y = 0; y < half_height; ++y) {
        memcpy(temp2.get(), upper, byte_pitch);
        memcpy(upper, temp.get(), byte_pitch);
        lower += byte_pitch;
        memcpy(temp.get(), lower, byte_pitch);
        memcpy(lower, temp2.get(), byte_pitch);
        upper += byte_pitch;
      }

      memcpy(upper, temp.get(), byte_pitch);
    }
  } else {
    auto* upper = (uint8_t*)py->GetData();
    auto* lower =
        (uint8_t*)(py->GetData() + (int64_t)half_height * py->GetPitch());
    for (int y = 0; y < half_height; ++y) {
      memcpy(temp.get(), upper, byte_pitch);
      memcpy(upper, lower, byte_pitch);
      memcpy(lower, temp.get(), byte_pitch);
      upper += byte_pitch;
      lower += byte_pitch;
    }
  }
}

void SwapV(Pyramid* py) { SwapUnswapV(py, false); }

void UnswapV(Pyramid* py) { SwapUnswapV(py, true); }

void Timer::Report(const char* name) {
  Output(1, "{}: {:.3f}s", name, Read());
};

}  // namespace multiblend::utils
