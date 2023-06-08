// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mb/multiblend.h"

#include <cmath>
#include <memory>
#include <optional>
#include <vector>

#include "mb/flex.h"
#include "mb/functions.h"
#include "mb/image.h"
#include "mb/linux_overrides.h"
#include "mb/logging.h"
#include "mb/pnger.h"

namespace multiblend {

namespace {

struct RecordState {
  int mc = 0;
  int prev_i = -1;
  int current_i = -1;
};

void Record(int tmp, int count, int x, RecordState& state,
            std::vector<io::Image>& images,
            std::optional<io::png::Pnger>& seam_map) {
  if (tmp == state.current_i) {
    state.mc += count;
    return;
  }

  int n_images = (int)images.size();

  if (state.mc > 0) {
    if (seam_map) {
      memset(&seam_map->line_[x - state.mc], state.current_i, state.mc);
    }
    for (int i = 0; i < n_images; ++i) {
      if (i == state.current_i) {
        images[i].masks_[0].Write32(0xc0000000 | state.mc);
      } else if (i == state.prev_i || state.prev_i == -1) {
        images[i].masks_[0].Write32(0x80000000 | state.mc);
      } else {
        images[i].masks_[0].IncrementLast32(state.mc);
      }
    }
  }

  state.prev_i = state.current_i;
  state.mc = count;
  state.current_i = tmp;
}

struct PyramidWithMasks : public Pyramid {
 public:
  using Pyramid::Pyramid;
  std::vector<utils::Flex> masks;
};
}  // namespace

Result Multiblend(std::vector<io::Image>& images, Options opts,
                  mt::ThreadpoolPtr threadpool) {
  utils::Timer timer;
  timer.Start();
  TimingResult timing;

  int n_images = (int)images.size();

  auto maskval = [&images](uint64_t x) {
    return (x & 0x7fffffffffffffff) | images[x & 0xffffffff].mask_state_;
  };

  /***********************************************************************
   * Open images to get prelimary info
   ***********************************************************************/
  std::size_t untrimmed_bytes = 0;

  for (int i = 0; i < n_images; ++i) {
    images[i].Open();
    untrimmed_bytes = std::max(untrimmed_bytes, images[i].untrimmed_bytes_);
  }

  /***********************************************************************
   * Check paramters, display warnings
   ***********************************************************************/
  for (int i = 1; i < n_images; ++i) {
    if (images[i].tiff_xres_ != images[0].tiff_xres_ ||
        images[i].tiff_yres_ != images[0].tiff_yres_) {
      utils::Output(0, "Warning: TIFF resolution mismatch ({} {}/{} {})",
                    images[0].tiff_xres_, images[0].tiff_yres_,
                    images[i].tiff_xres_, images[i].tiff_yres_);
    }
  }

  for (int i = 0; i < n_images; ++i) {
    if (opts.output_bpp == 0 && images[i].bpp_ == 16) {
      opts.output_bpp = 16;
    }
    if (images[i].bpp_ != images[0].bpp_) {
      utils::Throw(
          "Error: mixture of 8bpp and 16bpp images detected (not currently "
          "handled)");
    }
  }

  if (opts.output_bpp == 0) {
    opts.output_bpp = 8;
  } else if (opts.output_bpp == 16 &&
             opts.output_type == io::ImageType::MB_JPEG) {
    utils::Output(0, "Warning: 8bpp output forced by JPEG output");
    opts.output_bpp = 8;
  }

  /***********************************************************************
   * Allocate working space for reading/trimming/extraction
   ***********************************************************************/
  auto untrimmed_data = memory::AllocAligned<void>(untrimmed_bytes);

  /***********************************************************************
   * Read/trim/extract
   ***********************************************************************/
  for (int i = 0; i < n_images; ++i) {
    images[i].Read(untrimmed_data.get(), opts.gamma, threadpool);
  }

  /***********************************************************************
   * Clean up
   ***********************************************************************/
  untrimmed_data.reset();

  /***********************************************************************
   * Tighten
   ***********************************************************************/
  int min_xpos = 0x7fffffff;
  int min_ypos = 0x7fffffff;
  int width = 0;
  int height = 0;

  for (int i = 0; i < n_images; ++i) {
    min_xpos = std::min(min_xpos, images[i].xpos_);
    min_ypos = std::min(min_ypos, images[i].ypos_);
  }

  for (int i = 0; i < n_images; ++i) {
    images[i].xpos_ -= min_xpos;
    images[i].ypos_ -= min_ypos;
    width = std::max(width, images[i].xpos_ + images[i].width_);
    height = std::max(height, images[i].ypos_ + images[i].height_);
  }

  timing.images_time = timer.Read();

  /***********************************************************************
   * Determine number of levels
   ***********************************************************************/
  int blend_wh;
  int blend_levels;

  if (opts.fixed_levels == 0) {
    if (!opts.wideblend) {
      std::vector<int> widths;
      std::vector<int> heights;

      for (const auto& image : images) {
        widths.push_back(image.width_);
        heights.push_back(image.height_);
      }

      std::sort(widths.begin(), widths.end());
      std::sort(heights.begin(), heights.end());

      std::size_t halfway = (widths.size() - 1) >> 1;

      blend_wh =
          std::max((widths.size() & 1) != 0u
                       ? widths[halfway]
                       : (widths[halfway] + widths[halfway + 1] + 1) >> 1,
                   (heights.size() & 1) != 0u
                       ? heights[halfway]
                       : (heights[halfway] + heights[halfway + 1] + 1) >> 1);
    } else {
      blend_wh = (std::max)(width, height);
    }

    blend_levels = (int)std::floor(std::log2(blend_wh + 4.0f) - 1);
    if (opts.wideblend) {
      blend_levels++;
    }
  } else {
    blend_levels = opts.fixed_levels;
  }

  blend_levels += opts.add_levels;

  if (n_images == 1) {
    blend_levels = 0;
    utils::Output(1, "{} x {}, {} bpp", width, height, opts.output_bpp);
  } else {
    utils::Output(1, "{} x {}, {} levels, {} bpp", width, height, blend_levels,
                  opts.output_bpp);
  }

  /***********************************************************************
  ************************************************************************
  * Seaming
  ************************************************************************
  ***********************************************************************/
  timer.Start();

  utils::Output(1, "Seaming");
  switch ((static_cast<int>(!(opts.seamsave_filename == nullptr)) << 1) |
          static_cast<int>(!(opts.xor_filename == nullptr))) {
    case 1:
      utils::Output(1, " (saving XOR map)");
      break;
    case 2:
      utils::Output(1, " (saving seam map)");
      break;
    case 3:
      utils::Output(1, " (saving XOR and seam maps)");
      break;
  }

  int min_count;
  int xor_count;
  int xor_image;
  uint64_t utemp;
  int stop;

  uint64_t best;
  uint64_t a;
  uint64_t b;
  uint64_t c;
  uint64_t d;

  const uint64_t dt_max = 0x9000000000000000;
  uint64_t* prev_line = nullptr;
  uint64_t* this_line = nullptr;
  bool last_pixel = false;
  bool arbitrary_seam = false;

  auto seam_flex = utils::Flex(width, height);
  int max_queue = 0;

  /***********************************************************************
   * Backward distance transform
   ***********************************************************************/

  int n_threads = std::max(2, threadpool->GetNThreads());
  std::vector<std::vector<uint64_t>> thread_lines(n_threads);

  if (opts.seamload_filename == nullptr) {
    std::vector<std::vector<uint8_t>> thread_comp_lines(n_threads);

    for (int i = 0; i < n_threads; ++i) {
      thread_lines[i].resize(width);
      thread_comp_lines[i].resize(width);
    }

    // set all image masks to bottom right
    for (int i = 0; i < n_images; ++i) {
      images[i].tiff_mask_->End();
    }

    auto tasks = mt::MultiFuture<std::pair<uint8_t*, int>>{};
    for (int y = height - 1; y >= 0; --y) {
      int t = y % n_threads;
      this_line = thread_lines[t].data();
      uint8_t* comp = thread_comp_lines[t].data();

      // set initial image mask states
      for (int i = 0; i < n_images; ++i) {
        images[i].mask_state_ = 0x8000000000000000;
        if (y >= images[i].ypos_ && y < images[i].ypos_ + images[i].height_) {
          images[i].mask_count_ = width - (images[i].xpos_ + images[i].width_);
          images[i].mask_limit_ = images[i].xpos_;
        } else {
          images[i].mask_count_ = width;
          images[i].mask_limit_ = width;
        }
      }

      int x = width - 1;

      if (tasks.size() == n_threads) {
        tasks.wait();
        for (auto [comp_line, length] : tasks.get()) {
          seam_flex.Copy(comp_line, length);
          seam_flex.NextLine();
        }
        tasks = {};
      }

      while (x >= 0) {
        min_count = x + 1;
        xor_count = 0;

        // update image mask states
        for (int i = 0; i < n_images; ++i) {
          if (images[i].mask_count_ == 0) {
            if (x >= images[i].mask_limit_) {
              utemp = images[i].tiff_mask_->ReadBackwards32();
              images[i].mask_state_ = ((~utemp) << 32) & 0x8000000000000000;
              images[i].mask_count_ = utemp & 0x7fffffff;
            } else {
              images[i].mask_state_ = 0x8000000000000000;
              images[i].mask_count_ = min_count;
            }
          }

          if (images[i].mask_count_ < min_count) {
            min_count = images[i].mask_count_;
          }
          if (images[i].mask_state_ == 0u) {  // mask_state is inverted
            ++xor_count;
            xor_image = i;
          }
        }

        stop = x - min_count;

        if (xor_count == 1) {
          images[xor_image].seam_present_ = true;
          while (x > stop) {
            this_line[x--] = xor_image;
          }
        } else {
          if (y == height - 1) {   // bottom row
            if (x == width - 1) {  // first pixel(s)
              while (x > stop) {
                this_line[x--] = dt_max;  // max
              }
            } else {
              utemp = this_line[x + 1];
              utemp = maskval(utemp);
              while (x > stop) {
                utemp += 0x300000000;
                this_line[x--] = utemp;  // was min(temp, DT_MAX) but this is
                                         // unlikely to happen
              }
            }
          } else {                 // other rows
            if (x == width - 1) {  // first pixel(s)
              utemp = prev_line[x - 1] + 0x400000000;
              a = maskval(utemp);

              utemp = prev_line[x] + 0x300000000;
              b = maskval(utemp);

              d = a < b ? a : b;

              this_line[x--] = d;

              if (x == stop) {
                for (int i = 0; i < n_images; ++i) {
                  images[i].mask_count_ -= min_count;
                }
                continue;
              }

              c = b + 0x100000000;
              b = a - 0x100000000;
              d += 0x300000000;
            } else {
              utemp = prev_line[x] + 0x300000000;
              b = maskval(utemp);

              utemp = prev_line[x + 1] + 0x400000000;
              c = maskval(utemp);

              utemp = this_line[x + 1] + 0x300000000;
              d = maskval(utemp);
            }

            if (stop == -1) {
              stop = 0;
              last_pixel = true;
            }

            while (x > stop) {
              utemp = prev_line[x - 1] + 0x400000000;
              a = maskval(utemp);

              if (a < d) {
                d = a;
              }
              if (b < d) {
                d = b;
              }
              if (c < d) {
                d = c;
              }

              this_line[x--] = d;

              c = b + 0x100000000;
              b = a - 0x100000000;
              d += 0x300000000;
            }

            if (last_pixel) {
              // d is the new "best" to compare against
              if (b < d) {
                d = b;
              }
              if (c < d) {
                d = c;
              }

              this_line[x--] = d;

              last_pixel = false;
            }
          }
        }

        for (int i = 0; i < n_images; ++i) {
          images[i].mask_count_ -= min_count;
        }
      }

      if (y != 0) {
        tasks.push_back(threadpool->Queue([=] {
          int p = utils::CompressSeamLine(this_line, comp, width);
          if (p > width) {
            utils::Throw("bad p: {} at line {}", p, y);
          }
          return std::pair{comp, p};
        }));
      }

      prev_line = this_line;
    }  // end of row loop

    tasks.wait();
    for (auto [comp_line, length] : tasks.get()) {
      seam_flex.Copy(comp_line, length);
      seam_flex.NextLine();
    }

    for (int i = 0; i < n_images; ++i) {
      if (!images[i].seam_present_) {
        utils::Output(1, "Warning: {} is fully obscured by other images",
                      images[i].filename_);
      }
    }

  } else {  // if seamload_filename:
    for (int i = 0; i < n_images; ++i) {
      images[i].tiff_mask_->Start();
    }
  }

  // create top level masks
  for (int i = 0; i < n_images; ++i) {
    images[i].masks_.emplace_back(width, height);
  }

  std::optional<io::png::Pnger> xor_map;
  if (opts.xor_filename != nullptr) {
    xor_map = io::png::Pnger(opts.xor_filename, "XOR map", width, height,
                             io::png::ColorType::PALETTE);
  }

  std::optional<io::png::Pnger> seam_map;
  if (opts.seamsave_filename != nullptr) {
    seam_map = io::png::Pnger(opts.seamsave_filename, "Seam map", width, height,
                              io::png::ColorType::PALETTE);
  }

  /***********************************************************************
   * Forward distance transform
   ***********************************************************************/
  int current_count = 0;
  int64_t current_step;
  uint64_t dt_val;

  prev_line = thread_lines[1].data();

  uint64_t total_pixels = 0;
  uint64_t channel_totals[3] = {0};

  utils::Flex full_mask(width, height);
  utils::Flex xor_mask(width, height);

  bool alpha = false;

  for (int y = 0; y < height; ++y) {
    for (int i = 0; i < n_images; ++i) {
      images[i].mask_state_ = 0x8000000000000000;
      if (y >= images[i].ypos_ && y < images[i].ypos_ + images[i].height_) {
        images[i].mask_count_ = images[i].xpos_;
        images[i].mask_limit_ = images[i].xpos_ + images[i].width_;
      } else {
        images[i].mask_count_ = width;
        images[i].mask_limit_ = width;
      }
    }

    int x = 0;
    int best_temp;
    RecordState state;

    while (x < width) {
      min_count = width - x;
      xor_count = 0;

      for (int i = 0; i < n_images; ++i) {
        if (images[i].mask_count_ == 0) {
          if (x < images[i].mask_limit_) {
            utemp = images[i].tiff_mask_->ReadForwards32();
            images[i].mask_state_ = ((~utemp) << 32) & 0x8000000000000000;
            images[i].mask_count_ = utemp & 0x7fffffff;
          } else {
            images[i].mask_state_ = 0x8000000000000000;
            images[i].mask_count_ = min_count;
          }
        }

        if (images[i].mask_count_ < min_count) {
          min_count = images[i].mask_count_;
        }
        if (images[i].mask_state_ == 0u) {
          ++xor_count;
          xor_image = i;
        }
      }

      stop = x + min_count;

      if (xor_count == 0) {
        alpha = true;
      }
      full_mask.MaskWrite(min_count, xor_count != 0);
      xor_mask.MaskWrite(min_count, xor_count == 1);

      if (xor_count == 1) {
        if (xor_map) {
          memset(&xor_map->line_[x], xor_image, min_count);
        }

        std::size_t p =
            (y - images[xor_image].ypos_) * images[xor_image].width_ +
            (x - images[xor_image].xpos_);

        int total_count = min_count;
        total_pixels += total_count;
        if (opts.gamma) {
          switch (images[xor_image].bpp_) {
            case 8: {
              uint16_t v;
              while ((total_count--) != 0) {
                v = ((uint8_t*)images[xor_image].channels_[0].data_.get())[p];
                channel_totals[0] += static_cast<uint64_t>(v) * v;
                v = ((uint8_t*)images[xor_image].channels_[1].data_.get())[p];
                channel_totals[1] += static_cast<uint64_t>(v) * v;
                v = ((uint8_t*)images[xor_image].channels_[2].data_.get())[p];
                channel_totals[2] += static_cast<uint64_t>(v) * v;
                ++p;
              }
            } break;
            case 16: {
              uint32_t v;
              while ((total_count--) != 0) {
                v = ((uint16_t*)images[xor_image].channels_[0].data_.get())[p];
                channel_totals[0] += static_cast<uint64_t>(v) * v;
                v = ((uint16_t*)images[xor_image].channels_[1].data_.get())[p];
                channel_totals[1] += static_cast<uint64_t>(v) * v;
                v = ((uint16_t*)images[xor_image].channels_[2].data_.get())[p];
                channel_totals[2] += static_cast<uint64_t>(v) * v;
                ++p;
              }
            } break;
          }
        } else {
          switch (images[xor_image].bpp_) {
            case 8: {
              while ((total_count--) != 0) {
                channel_totals[0] +=
                    ((uint8_t*)images[xor_image].channels_[0].data_.get())[p];
                channel_totals[1] +=
                    ((uint8_t*)images[xor_image].channels_[1].data_.get())[p];
                channel_totals[2] +=
                    ((uint8_t*)images[xor_image].channels_[2].data_.get())[p];
                ++p;
              }
            } break;
            case 16: {
              while ((total_count--) != 0) {
                channel_totals[0] +=
                    ((uint16_t*)images[xor_image].channels_[0].data_.get())[p];
                channel_totals[1] +=
                    ((uint16_t*)images[xor_image].channels_[1].data_.get())[p];
                channel_totals[2] +=
                    ((uint16_t*)images[xor_image].channels_[2].data_.get())[p];
                ++p;
              }
            } break;
          }
        }

        if (opts.seamload_filename == nullptr) {
          Record(xor_image, min_count, x, state, images, seam_map);
          while (x < stop) {
            this_line[x++] = xor_image;
          }
        } else {
          x = stop;
        }

        best = xor_image;
      } else {
        if (xor_map) {
          memset(&xor_map->line_[x], 0xff, min_count);
        }

        if (opts.seamload_filename == nullptr) {
          if (y == 0) {
            // top row
            while (x < stop) {
              best = this_line[x];

              if (x > 0) {
                utemp = this_line[x - 1] + 0x300000000;
                d = maskval(utemp);

                if (d < best) {
                  best = d;
                }
              }

              if (((best & 0x8000000000000000) != 0u) && (xor_count != 0)) {
                arbitrary_seam = true;
                for (int i = 0; i < n_images; ++i) {
                  if (images[i].mask_state_ == 0u) {
                    best = 0x8000000000000000 | i;
                    if (!opts.reverse) {
                      break;
                    }
                  }
                }
              }

              best_temp = best & 0xffffffff;
              Record(best_temp, 1, x, state, images, seam_map);
              this_line[x++] = best;
            }
          } else {
            // other rows
            if (x == 0) {
              utils::ReadSeamDT(seam_flex, current_count, current_step, dt_val);
              best = dt_val;

              utemp = *prev_line + 0x300000000;
              b = maskval(utemp);
              if (b < best) {
                best = b;
              }

              utemp = prev_line[1] + 0x400000000;
              c = maskval(utemp);
              if (c < best) {
                best = c;
              }

              if (((best & 0x8000000000000000) != 0u) && (xor_count != 0)) {
                arbitrary_seam = true;
                for (int i = 0; i < n_images; ++i) {
                  if (images[i].mask_state_ == 0u) {
                    best = 0x8000000000000000 | i;
                    if (!opts.reverse) {
                      break;
                    }
                  }
                }
              }

              best_temp = best & 0xffffffff;
              Record(best_temp, 1, x, state, images, seam_map);
              this_line[x++] = best;

              if (x == stop) {
                for (int i = 0; i < n_images; ++i) {
                  images[i].mask_count_ -= min_count;
                }
                continue;
              }

              a = b + 0x100000000;
              b = c - 0x100000000;
            } else {
              utemp = prev_line[x - 1] + 0x400000000;
              a = maskval(utemp);

              utemp = prev_line[x] + 0x300000000;
              b = maskval(utemp);
            }

            utemp = best + 0x300000000;
            d = maskval(utemp);

            if (stop == width) {
              stop--;
              last_pixel = true;
            }

            while (x < stop) {
              utemp = prev_line[x + 1] + 0x400000000;
              c = maskval(utemp);

              utils::ReadSeamDT(seam_flex, current_count, current_step, dt_val);
              best = dt_val;

              if (a < best) {
                best = a;
              }
              if (b < best) {
                best = b;
              }
              if (c < best) {
                best = c;
              }
              if (d < best) {
                best = d;
              }

              if (((best & 0x8000000000000000) != 0u) && (xor_count != 0)) {
                arbitrary_seam = true;
                for (int i = 0; i < n_images; ++i) {
                  if (images[i].mask_state_ == 0u) {
                    best = 0x8000000000000000 | i;
                    if (!opts.reverse) {
                      break;
                    }
                  }
                }
              }

              best_temp = best & 0xffffffff;
              Record(best_temp, 1, x, state, images, seam_map);
              this_line[x++] = best;  // best;

              a = b + 0x100000000;
              b = c - 0x100000000;
              d = best + 0x300000000;
            }

            if (last_pixel) {
              utils::ReadSeamDT(seam_flex, current_count, current_step, dt_val);
              best = dt_val;

              if (a < best) {
                best = a;
              }
              if (b < best) {
                best = b;
              }
              if (d < best) {
                best = d;
              }

              if (((best & 0x8000000000000000) != 0u) && (xor_count != 0)) {
                arbitrary_seam = true;
                for (int i = 0; i < n_images; ++i) {
                  if (images[i].mask_state_ == 0u) {
                    best = 0x8000000000000000 | i;
                    if (!opts.reverse) {
                      break;
                    }
                  }
                }
              }

              best_temp = best & 0xffffffff;
              Record(best_temp, 1, x, state, images, seam_map);
              this_line[x++] = best;  // best;

              last_pixel = false;
            }
          }
        } else {  // if (seamload_filename)...
          x = stop;
        }
      }

      for (int i = 0; i < n_images; ++i) {
        images[i].mask_count_ -= min_count;
      }
    }

    if (opts.seamload_filename == nullptr) {
      Record(-1, 0, x, state, images, seam_map);

      for (int i = 0; i < n_images; ++i) {
        images[i].masks_[0].NextLine();
      }
    }

    full_mask.NextLine();
    xor_mask.NextLine();

    if (xor_map) {
      xor_map->Write();
    }
    if (seam_map) {
      seam_map->Write();
    }

    std::swap(this_line, prev_line);
  }

  thread_lines.resize(0);
  xor_map.reset();
  seam_map.reset();

  if (!alpha || opts.output_type == io::ImageType::MB_JPEG) {
    opts.no_mask = true;
  }

/***********************************************************************
 * Seam load
 ***********************************************************************/
#ifdef MULTIBLEND_WITH_PNG
  if (opts.seamload_filename != nullptr) {
    int png_depth;
    int png_colour;
    png_uint_32 png_width;
    png_uint_32 png_height;
    uint8_t sig[8];

    FILE* tmp_file;
    fopen_s(&tmp_file, opts.seamload_filename, "rb");
    if (tmp_file == nullptr) {
      utils::Throw("Error: Couldn't open seam file");
    }
    auto output_file = std::unique_ptr<FILE, io::FileDeleter>{tmp_file};

    std::size_t r =
        fread(sig, 1, 8,
              output_file.get());  // assignment suppresses g++ -Ofast warning
    if (!png_check_sig(sig, 8)) {
      utils::Throw("Error: Bad PNG signature");
    }

    auto png_ptr = std::unique_ptr<png_struct, io::png::PngReadStructDeleter>{
        png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr,
                               nullptr),
        io::png::PngReadStructDeleter{}};

    if (png_ptr == nullptr) {
      utils::Throw("Error: Seam PNG problem");
    }

    auto info_ptr = std::unique_ptr<png_info, io::png::PngInfoStructDeleter>{
        png_create_info_struct(png_ptr.get()),
        io::png::PngInfoStructDeleter{png_ptr.get()}};

    if (info_ptr == nullptr) {
      utils::Throw("Error: Seam PNG problem");
    }

    png_init_io(png_ptr.get(), output_file.get());
    png_set_sig_bytes(png_ptr.get(), 8);
    png_read_info(png_ptr.get(), info_ptr.get());
    png_get_IHDR(png_ptr.get(), info_ptr.get(), &png_width, &png_height,
                 &png_depth, &png_colour, nullptr, nullptr, nullptr);

    if (png_width != width || png_height != height) {
      utils::Throw("Error: Seam PNG dimensions don't match workspace");
    }
    if (png_depth != 8 || png_colour != PNG_COLOR_TYPE_PALETTE) {
      utils::Throw("Error: Incorrect seam PNG format");
    }

    auto png_line = std::make_unique<png_byte[]>(width);
    std::optional<io::png::Pnger> empty;

    for (int y = 0; y < height; ++y) {
      png_read_row(png_ptr.get(), png_line.get(), nullptr);

      int x = 0;
      RecordState state;

      for (x = 0; x < width; ++x) {
        if (png_line[x] > n_images) {
          utils::Throw("Error: Bad pixel found in seam file: {},{}", x, y);
        }
        Record(png_line[x], 1, x, state, images, empty);
      }

      Record(-1, 0, x, state, images, empty);

      for (int i = 0; i < n_images; ++i) {
        images[i].masks_[0].NextLine();
      }
    }
  }
#endif

  timing.seam_time = timer.Read();

  /***********************************************************************
   * No output?
   ***********************************************************************/
  std::array<memory::AlignedAllocPtr<void>, 3> output_channels = {
      nullptr, nullptr, nullptr};

  if (opts.output_type != io::ImageType::MB_NONE) {
    /***********************************************************************
     * Shrink masks
     ***********************************************************************/
    utils::Output(1, "Shrinking masks...");

    timer.Start();
    {
      auto tasks = mt::MultiFuture<void>{};
      for (int i = 0; i < n_images; ++i) {
        tasks.push_back(threadpool->Queue(
            [=, &images] { ShrinkMasks(images[i].masks_, blend_levels); }));
      }
      tasks.wait();
      tasks.get();
    }

    timing.shrink_mask_time = timer.Read();

    /***********************************************************************
     * Create shared input pyramids
     ***********************************************************************/
    // wrapping
    std::vector<PyramidWithMasks> wrap_pyramids;
    int wrap_levels_h = 0;
    int wrap_levels_v = 0;

    if ((opts.wrap & 1) != 0) {
      wrap_levels_h = (int)std::floor(std::log2((width >> 1) + 4.0f) - 1);
      wrap_pyramids.emplace_back(width >> 1, height, wrap_levels_h, 0, 0,
                                 threadpool);
      wrap_pyramids.emplace_back((width + 1) >> 1, height, wrap_levels_h,
                                 width >> 1, 0, threadpool);
    }

    if ((opts.wrap & 2) != 0) {
      wrap_levels_v = (int)std::floor(std::log2((height >> 1) + 4.0f) - 1);
      wrap_pyramids.emplace_back(width, height >> 1, wrap_levels_v, 0, 0,
                                 threadpool);
      wrap_pyramids.emplace_back(width, (height + 1) >> 1, wrap_levels_v, 0,
                                 height >> 1, threadpool);
    }

    // masks
    {
      auto tasks = mt::MultiFuture<void>{};
      for (auto& pyr : wrap_pyramids) {
        tasks.push_back(threadpool->Queue([=, &pyr] {
          pyr.masks.emplace_back(width, height);
          for (int y = 0; y < height; ++y) {
            if (y < pyr.GetY() || y >= pyr.GetY() + pyr.GetHeight()) {
              pyr.masks[0].Write32(0x80000000 | width);
            } else {
              if (pyr.GetX() != 0) {
                pyr.masks[0].Write32(0x80000000 | pyr.GetX());
                pyr.masks[0].Write32(0xc0000000 | pyr.GetWidth());
              } else {
                pyr.masks[0].Write32(0xc0000000 | pyr.GetWidth());
                if (pyr.GetWidth() != width) {
                  pyr.masks[0].Write32(0x80000000 | (width - pyr.GetWidth()));
                }
              }
            }
            pyr.masks[0].NextLine();
          }

          ShrinkMasks(pyr.masks,
                      pyr.GetWidth() == width ? wrap_levels_v : wrap_levels_h);
        }));
      }
      tasks.wait();
      tasks.get();
    }
    // end wrapping

    int total_levels =
        std::max({blend_levels, wrap_levels_h, wrap_levels_v, 1});

    for (int i = 0; i < n_images; ++i) {
      images[i].pyramid_ =
          Pyramid(images[i].width_, images[i].height_, blend_levels,
                  images[i].xpos_, images[i].ypos_, threadpool);
    }

    for (int level = total_levels - 1; level >= 0; --level) {
      std::size_t max_bytes = 0;

      if (level < blend_levels) {
        for (const auto& image : images) {
          max_bytes =
              std::max(max_bytes, image.pyramid_->GetLevel(level).bytes);
        }
      }

      for (auto& pyr : wrap_pyramids) {
        if (level < pyr.GetNLevels()) {
          max_bytes = std::max(max_bytes, pyr.GetLevel(level).bytes);
        }
      }

      auto temp =
          std::shared_ptr<float>(memory::AllocAligned<float>(max_bytes));

      if (level < blend_levels) {
        for (auto& image : images) {
          image.pyramid_->GetLevel(level).data = temp;
        }
      }

      for (auto& pyr : wrap_pyramids) {
        if (level < pyr.GetNLevels()) {
          pyr.GetLevel(level).data = temp;
        }
      }
    }

    /***********************************************************************
     * Create output pyramid
     ***********************************************************************/
    auto output_pyramid =
        Pyramid{width, height, total_levels, 0, 0, threadpool};

    for (int level = total_levels - 1; level >= 0; --level) {
      output_pyramid.GetLevel(level).data =
          memory::AllocAligned<float>(output_pyramid.GetLevel(level).bytes);
    }

    /***********************************************************************
     * Blend
     ***********************************************************************/
    if (n_images == 1) {
      if (opts.wrap != 0) {
        utils::Output(1, "Wrapping...");
      } else {
        utils::Output(1, "Processing...");
      }
    } else {
      if (opts.wrap != 0) {
        utils::Output(1, "Blending/wrapping...");
      } else {
        utils::Output(1, "Blending...");
      }
    }

    for (int c = 0; c < 3; ++c) {
      if (n_images > 1) {
        for (int i = 0; i < n_images; ++i) {
          timer.Start();

          images[i].pyramid_->Copy((uint8_t*)images[i].channels_[c].data_.get(),
                                   1, images[i].width_, opts.gamma,
                                   images[i].bpp_);
          if (opts.output_bpp != images[i].bpp_) {
            images[i].pyramid_->Multiply(
                0, opts.gamma ? (opts.output_bpp == 8 ? 1.0f / 66049 : 66049)
                              : (opts.output_bpp == 8 ? 1.0f / 257 : 257));
          }

          images[i].channels_[c].data_.reset();

          timing.copy_time += timer.Read();

          timer.Start();
          images[i].pyramid_->Shrink();
          timing.shrink_time += timer.Read();

          timer.Start();
          images[i].pyramid_->Laplace();
          timing.laplace_time += timer.Read();

          // blend into output pyramid...

          timer.Start();

          for (int level = 0; level < blend_levels; ++level) {
            const auto& in_level = images[i].pyramid_->GetLevel(level);
            auto& out_level = output_pyramid.GetLevel(level);

            int x_offset = (in_level.x - out_level.x) >> level;
            int y_offset = (in_level.y - out_level.y) >> level;

            auto tasks = mt::MultiFuture<void>{};
            for (int b = 0; b < (int)out_level.bands.size() - 1; ++b) {
              int sy = out_level.bands[b];
              int ey = out_level.bands[b + 1];

              tasks.push_back(
                  threadpool->Queue([=, &images, &in_level, &out_level] {
                    for (int y = sy; y < ey; ++y) {
                      int in_line = y - y_offset;
                      if (in_line < 0) {
                        in_line = 0;
                      } else if (in_line > in_level.height - 1) {
                        in_line = in_level.height - 1;
                      }
                      float* input_p = in_level.data.get() +
                                       (std::size_t)in_line * in_level.pitch;
                      float* output_p = out_level.data.get() +
                                        (std::size_t)y * out_level.pitch;

                      utils::CompositeLine(input_p, output_p, i, x_offset,
                                           in_level.width, out_level.width,
                                           out_level.pitch,
                                           images[i].masks_[level].data_.get(),
                                           images[i].masks_[level].rows_[y]);
                    }
                  }));
            }
            tasks.wait();
            tasks.get();
          }

          timing.blend_time += timer.Read();
        }

        timer.Start();
        output_pyramid.Collapse(blend_levels);
        timing.collapse_time += timer.Read();
      } else {
        timer.Start();

        output_pyramid.Copy((uint8_t*)images[0].channels_[c].data_.get(), 1,
                            images[0].width_, opts.gamma, images[0].bpp_);
        if (opts.output_bpp != images[0].bpp_) {
          output_pyramid.Multiply(
              0, opts.gamma ? (opts.output_bpp == 8 ? 1.0f / 66049 : 66049)
                            : (opts.output_bpp == 8 ? 1.0f / 257 : 257));
        }

        images[0].channels_[c].data_.reset();

        timing.copy_time += timer.Read();
      }

      /***********************************************************************
       * Wrapping
       ***********************************************************************/
      if (opts.wrap != 0) {
        timer.Start();

        int p = 0;

        for (int w = 1; w <= 2; ++w) {
          if ((opts.wrap & w) != 0) {
            if (w == 1) {
              utils::SwapH(&output_pyramid);
            } else {
              utils::SwapV(&output_pyramid);
            }

            int wrap_levels = (w == 1) ? wrap_levels_h : wrap_levels_v;
            for (int wp = 0; wp < 2; ++wp) {
              wrap_pyramids[p].Copy(
                  (uint8_t*)(output_pyramid.GetData() +
                             wrap_pyramids[p].GetX() +
                             wrap_pyramids[p].GetY() *
                                 (int64_t)output_pyramid.GetPitch()),
                  1, output_pyramid.GetPitch(), false, 32);
              wrap_pyramids[p].Shrink();
              wrap_pyramids[p].Laplace();

              for (int level = 0; level < wrap_levels; ++level) {
                const auto& in_level = wrap_pyramids[p].GetLevel(level);
                auto& out_level = output_pyramid.GetLevel(level);

                int x_offset = (in_level.x - out_level.x) >> level;
                int y_offset = (in_level.y - out_level.y) >> level;

                auto tasks = mt::MultiFuture<void>{};
                for (int b = 0; b < (int)out_level.bands.size() - 1; ++b) {
                  int sy = out_level.bands[b];
                  int ey = out_level.bands[b + 1];

                  tasks.push_back(threadpool->Queue([=, &in_level, &out_level,
                                                     &wrap_pyramids] {
                    for (int y = sy; y < ey; ++y) {
                      int in_line = y - y_offset;
                      if (in_line < 0) {
                        in_line = 0;
                      } else if (in_line > in_level.height - 1) {
                        in_line = in_level.height - 1;
                      }
                      float* input_p = in_level.data.get() +
                                       (std::size_t)in_line * in_level.pitch;
                      float* output_p = out_level.data.get() +
                                        (std::size_t)y * out_level.pitch;

                      utils::CompositeLine(
                          input_p, output_p, wp + static_cast<int>(level == 0),
                          x_offset, in_level.width, out_level.width,
                          out_level.pitch,
                          wrap_pyramids[p].masks[level].data_.get(),
                          wrap_pyramids[p].masks[level].rows_[y]);
                    }
                  }));
                }
                tasks.wait();
                tasks.get();
              }
              ++p;
            }

            output_pyramid.Collapse(wrap_levels);

            if (w == 1) {
              utils::UnswapH(&output_pyramid);
            } else {
              utils::UnswapV(&output_pyramid);
            }
          }  // if (wrap & w)
        }    // w loop

        timing.wrap_time += timer.Read();
      }  // if (wrap)
         // end wrapping

      /***********************************************************************
       * Offset correction
       ***********************************************************************/
      if (total_pixels != 0u) {
        double channel_total = 0;  // must be a double
        float* data = output_pyramid.GetData();
        xor_mask.Start();

        for (int y = 0; y < height; ++y) {
          int x = 0;
          while (x < width) {
            uint32_t v = xor_mask.ReadForwards32();
            if ((v & 0x80000000) != 0u) {
              v = x + v & 0x7fffffff;
              while (x < (int)v) {
                channel_total += data[x++];
              }
            } else {
              x += v;
            }
          }

          data += output_pyramid.GetPitch();
        }

        float avg = (float)channel_totals[c] / total_pixels;
        if (opts.output_bpp != images[0].bpp_) {
          switch (opts.output_bpp) {
            case 8:
              avg /= 256;
              break;
            case 16:
              avg *= 256;
              break;
          }
        }
        float output_avg = (float)channel_total / total_pixels;
        output_pyramid.Add(avg - output_avg, 1);
      }

      /***********************************************************************
       * Output
       ***********************************************************************/
      timer.Start();

      output_channels[c] = memory::AllocAligned<void>(
          ((std::size_t)width * height) << (opts.output_bpp >> 4));

      switch (opts.output_bpp) {
        case 8:
          output_pyramid.Out((uint8_t*)output_channels[c].get(), width,
                             opts.gamma, opts.dither, true);
          break;
        case 16:
          output_pyramid.Out((uint16_t*)output_channels[c].get(), width,
                             opts.gamma, opts.dither, true);
          break;
      }

      timing.out_time += timer.Read();
    }
  }

  full_mask.Shrink();

  return Result{
      .output_bpp = opts.output_bpp,
      .width = width,
      .height = height,
      .no_mask = opts.no_mask,

      .output_channels = std::move(output_channels),
      .min_xpos = min_xpos,
      .min_ypos = min_ypos,
      .full_mask = std::move(full_mask),

      .timing = timing,
  };
}

}  // namespace multiblend
