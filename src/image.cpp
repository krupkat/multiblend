#include "mb/image.h"

#include <cmath>
#include <cstdio>
#include <memory>
#include <optional>
#include <string>

#include <spdlog/fmt/fmt.h>

#include "mb/functions.h"
#ifdef MULTIBLEND_WITH_JPEG
#include "mb/jpeg.h"
#endif
#include "mb/linux_overrides.h"
#include "mb/logging.h"
#include "mb/mapalloc.h"
#include "mb/pnger.h"
#include "mb/pyramid.h"
#include "mb/threadpool.h"
#ifdef MULTIBLEND_WITH_TIFF
#include "mb/tiff.h"
#endif

namespace multiblend::io {

int hist_red[256];
int hist_grn[256];
int hist_blu[256];

Image::Image(char* filename) : filename_(filename) {}

Image::Image(InMemoryImage image)
    : image_(std::move(image)), filename_("in memory image") {}

/***********************************************************************
************************************************************************
** Open
************************************************************************
***********************************************************************/
void Image::Open() {
  float tiff_xpos;
  float tiff_ypos;
  uint16_t compression;

  if (image_) {
    type_ = ImageType::MB_IN_MEMORY;
  } else {
    const char* ext = strrchr(filename_.c_str(), '.');
    if (ext == nullptr) {
      utils::Throw("Could not identify file extension: {}", filename_);
    }
    ++ext;

    if ((_stricmp(ext, "tif") == 0) || (_stricmp(ext, "tiff") == 0)) {
      type_ = ImageType::MB_TIFF;
    } else if ((_stricmp(ext, "jpg") == 0) || (_stricmp(ext, "jpeg") == 0)) {
      type_ = ImageType::MB_JPEG;
    } else if (_stricmp(ext, "png") == 0) {
      type_ = ImageType::MB_PNG;
    } else {
      utils::Throw("Unknown file extension: {}", filename_);
    }
  }

  switch (type_) {
    case ImageType::MB_TIFF: {
#ifdef MULTIBLEND_WITH_TIFF
      tiff_ = {TIFFOpen(filename_.c_str(), "r"), tiff::CloseDeleter{}};
      if (tiff_ == nullptr) {
        utils::Throw("Could not open {}", filename_);
      }

      if (TIFFGetField(tiff_.get(), TIFFTAG_XPOSITION, &tiff_xpos) == 0) {
        tiff_xpos = -1;
      }
      if (TIFFGetField(tiff_.get(), TIFFTAG_YPOSITION, &tiff_ypos) == 0) {
        tiff_ypos = -1;
      }
      if (TIFFGetField(tiff_.get(), TIFFTAG_XRESOLUTION, &tiff_xres_) == 0) {
        tiff_xres_ = -1;
      }
      if (TIFFGetField(tiff_.get(), TIFFTAG_YRESOLUTION, &tiff_yres_) == 0) {
        tiff_yres_ = -1;
      }
      TIFFGetField(tiff_.get(), TIFFTAG_IMAGEWIDTH, &tiff_width_);
      TIFFGetField(tiff_.get(), TIFFTAG_IMAGELENGTH, &tiff_height_);
      TIFFGetField(tiff_.get(), TIFFTAG_BITSPERSAMPLE, &bpp_);
      TIFFGetField(tiff_.get(), TIFFTAG_SAMPLESPERPIXEL, &spp_);
      TIFFGetField(tiff_.get(), TIFFTAG_ROWSPERSTRIP, &rows_per_strip_);
      TIFFGetField(tiff_.get(), TIFFTAG_COMPRESSION, &compression);

      if (bpp_ != 8 && bpp_ != 16) {
        utils::Output(1, "Invalid bpp {} ({})", bpp_, filename_);
        utils::Output(1, "{}, {}", tiff_width_, tiff_height_);
        TIFFGetField(tiff_.get(), TIFFTAG_BITSPERSAMPLE, &bpp_);
        if (bpp_ != 8 && bpp_ != 16) {
          utils::Throw("Invalid bpp {} ({})", bpp_, filename_);
        }
      }
      //   if (spp != 4) die("Images must be RGBA (%s)",
      // filename);

      geotiff_.set = false;

      if (tiff_xpos == -1 && tiff_ypos == -1) {
        // try to read geotiff tags
        if (tiff::geotiff_read(tiff_.get(), &geotiff_) != 0) {
          xpos_ = (int)(geotiff_.XGeoRef / geotiff_.XCellRes);
          ypos_ = (int)(geotiff_.YGeoRef / geotiff_.YCellRes);
        } else {
          xpos_ = 0;
          ypos_ = 0;
        }
      } else {
        if (tiff_xpos != -1 && tiff_xres_ > 0) {
          xpos_ = (int)std::lroundf(tiff_xpos * tiff_xres_);
        }
        if (tiff_ypos != -1 && tiff_yres_ > 0) {
          ypos_ = (int)std::lroundf(tiff_ypos * tiff_yres_);
        }
      }

      first_strip_ = 0;
      end_strip_ = TIFFNumberOfStrips(tiff_.get());

      int s;
      tmsize_t temp;

      if (tiff_xpos <= 0 && tiff_ypos <= 0 && compression != 1 && spp_ == 4) {
        tmsize_t min_stripsize = 0xffffffff;
        int min_stripcount = 0;
        for (s = 0; s < (int)TIFFNumberOfStrips(tiff_.get()) - 1; ++s) {
          temp = TIFFRawStripSize(tiff_.get(), s);
          if (temp < min_stripsize) {
            min_stripsize = temp;
            min_stripcount = 1;
          } else if (temp == min_stripsize) {
            min_stripcount++;
          }
        }

        if (min_stripcount > 2) {
          first_strip_ = -1;
          for (s = 0; s < (int)TIFFNumberOfStrips(tiff_.get()); ++s) {
            temp = TIFFRawStripSize(tiff_.get(), s);
            if (temp != min_stripsize) {
              if (first_strip_ == -1) {
                first_strip_ = s;
              }
              end_strip_ = s + 1;
            }
          }
          if (first_strip_ == -1) {
            first_strip_ = 0;
          }
        }
      }

      if ((first_strip_ != 0) ||
          end_strip_ != TIFFNumberOfStrips(
                            tiff_.get())) {  // double check that min strips are
                                             // (probably) transparent
        auto buf = std::unique_ptr<void, tiff::FreeDeleter>{
            _TIFFmalloc(TIFFScanlineSize(tiff_.get())), tiff::FreeDeleter{}};
        if (first_strip_ != 0) {
          TIFFReadScanline(tiff_.get(), buf.get(), 0);
        } else {
          TIFFReadScanline(tiff_.get(), buf.get(), tiff_u_height_ - 1);
        }
        bool trans;
        switch (bpp_) {
          case 8:
            trans = (((uint32_t*)buf.get())[0] & 0xff000000) == 0u;
            break;
          case 16:
            trans = (((uint64_t*)buf.get())[0] & 0xffff000000000000) == 0u;
            break;
        }
        if (!trans) {
          first_strip_ = 0;
          end_strip_ = TIFFNumberOfStrips(tiff_.get());
        }
      }

      ypos_ += first_strip_ * rows_per_strip_;
      int rows_missing =
          TIFFNumberOfStrips(tiff_.get()) * rows_per_strip_ - tiff_height_;

      tiff_u_height_ = (end_strip_ - first_strip_) * rows_per_strip_;
      if (end_strip_ == TIFFNumberOfStrips(tiff_.get())) {
        tiff_u_height_ -= rows_missing;
      }
#else
      throw(std::runtime_error("TIFF support not compiled in"));
#endif
    } break;
    case ImageType::MB_JPEG: {
#ifdef MULTIBLEND_WITH_JPEG
      FILE* tmp_file = nullptr;
      fopen_s(&tmp_file, filename_.c_str(), "rb");
      if (tmp_file == nullptr) {
        utils::Throw("Could not open {}", filename_);
      }
      file_ = {tmp_file, FileDeleter{}};

      cinfo_ = {new jpeg_decompress_struct{}, jpeg::DecompressDeleter{}};
      jerr_ = std::make_unique<jpeg_error_mgr>();

      cinfo_->err = jpeg_std_error(jerr_.get());
      jpeg_create_decompress(cinfo_.get());
      jpeg_stdio_src(cinfo_.get(), file_.get());
      jpeg_read_header(cinfo_.get(), TRUE);
      jpeg_start_decompress(cinfo_.get());

      if ((cinfo_->output_width == 0u) || (cinfo_->output_height == 0u)) {
        utils::Throw("Unknown JPEG format ({})", filename_);
      }

      if (cinfo_->out_color_components != 3) {
        utils::Throw("Unknown JPEG format ({})", filename_);
      }

      tiff_width_ = cinfo_->output_width;
      tiff_height_ = tiff_u_height_ = cinfo_->output_height;

      bpp_ = 8;
      spp_ = 3;

      xpos_ = ypos_ = 0;
      tiff_xpos = tiff_ypos = 0;
      tiff_xres_ = tiff_yres_ = 90;
#else
      throw(std::runtime_error("JPEG support not compiled in"));
#endif
    } break;
    case ImageType::MB_PNG: {
#ifdef MULTIBLEND_WITH_PNG
      FILE* tmp_file = nullptr;
      fopen_s(&tmp_file, filename_.c_str(), "rb");
      if (tmp_file == nullptr) {
        utils::Throw("Could not open {}", filename_);
      }
      file_ = {tmp_file, FileDeleter{}};

      uint8_t sig[8];
      std::size_t r = fread(
          sig, 1, 8, file_.get());  // assignment suppresses g++ -Ofast warning
      if (!png_check_sig(sig, 8)) {
        utils::Throw("Bad PNG signature ({})", filename_);
      }

      png_ptr_ = {png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr,
                                         nullptr, nullptr),
                  png::PngReadStructDeleter{}};
      if (png_ptr_ == nullptr) {
        utils::Throw("Error: libpng problem");
      }

      auto info_ptr = std::unique_ptr<png_info, png::PngInfoStructDeleter>{
          png_create_info_struct(png_ptr_.get()),
          png::PngInfoStructDeleter{png_ptr_.get()}};

      if (info_ptr == nullptr) {
        utils::Throw("Error: libpng problem");
      }

      png_init_io(png_ptr_.get(), file_.get());
      png_set_sig_bytes(png_ptr_.get(), 8);
      png_read_info(png_ptr_.get(), info_ptr.get());

      int png_colour;
      uint32_t png_width;
      uint32_t png_height;
      int _bpp;
      png_get_IHDR(png_ptr_.get(), info_ptr.get(), &png_width, &png_height,
                   &_bpp, &png_colour, nullptr, nullptr, nullptr);
      bpp_ = _bpp;
      tiff_width_ = png_width;
      tiff_u_height_ = tiff_height_ = png_height;

      switch (png_colour) {
        case PNG_COLOR_TYPE_RGB:
          spp_ = 3;
          break;
        case PNG_COLOR_TYPE_RGBA:
          spp_ = 4;
          break;
        default:
          utils::Throw("Bad PNG colour type ({})", filename_);
      }

      if (bpp_ != 8 && bpp_ != 16) {
        utils::Throw("Bad bit depth ({})", filename_);
      }

      xpos_ = ypos_ = 0;
      tiff_xpos = tiff_ypos = 0;
      tiff_xres_ = tiff_yres_ = 90;
#else
      throw(std::runtime_error("PNG support not compiled in"));
#endif
    } break;
    case ImageType::MB_IN_MEMORY: {
      tiff_width_ = image_->tiff_width;
      tiff_height_ = tiff_u_height_ = image_->tiff_height;
      bpp_ = image_->bpp;
      spp_ = image_->spp;
      xpos_add_ = image_->xpos_add;
      ypos_add_ = image_->ypos_add;
      xpos_ = ypos_ = 0;
      tiff_xres_ = tiff_yres_ = 90;
    } break;
  }

  xpos_ += xpos_add_;
  ypos_ += ypos_add_;

  std::size_t untrimmed_pixels = (std::size_t)tiff_u_height_ * tiff_width_;
  untrimmed_bytes_ = (untrimmed_pixels * spp_) << (bpp_ >> 4);
}

/***********************************************************************
************************************************************************
** Read
************************************************************************
***********************************************************************/
void Image::Read(void* data, bool gamma) {
  utils::Output(1, "Processing {}...", filename_);

  switch (type_) {
    case ImageType::MB_TIFF: {
#ifdef MULTIBLEND_WITH_TIFF
      char* pointer = (char*)data;

      for (int s = first_strip_; s < end_strip_; s++) {
        auto strip_size = TIFFReadEncodedStrip(tiff_.get(), s, pointer, -1);
        pointer += strip_size;
      }
#else
      throw(std::runtime_error("TIFF support not compiled in"));
#endif
    } break;
    case ImageType::MB_JPEG: {
#ifdef MULTIBLEND_WITH_JPEG
      auto* pointer = (uint8_t*)data;

      while (cinfo_->output_scanline < cinfo_->output_height) {
        jpeg_read_scanlines(cinfo_.get(), &pointer, 1);
        pointer += (tiff_width_ * spp_) << (bpp_ >> 4);
      }
#else
      throw(std::runtime_error("JPEG support not compiled in"));
#endif
    } break;
    case ImageType::MB_PNG: {
#ifdef MULTIBLEND_WITH_PNG
      auto* pointer = (uint8_t*)data;

      for (int y = 0; y < tiff_height_; ++y) {
        png_read_row(png_ptr_.get(), pointer, nullptr);
        pointer += (tiff_width_ * spp_) << (bpp_ >> 4);
      }
#else
      throw(std::runtime_error("PNG support not compiled in"));
#endif
    } break;
    case ImageType::MB_IN_MEMORY: {
      auto size_bytes = image_->data.size() << (bpp_ >> 4);
      memcpy(data, image_->data.data(), size_bytes);
    } break;
  }

  /***********************************************************************
   * Trim
   ***********************************************************************/
  int x;
  int y;
  int top;
  int left;
  int bottom;
  int right;

  if (spp_ == 4) {
    switch (bpp_) {
      case 8: {
        auto* p = (uint32_t*)data;
        p--;

        for (y = 0; y < tiff_u_height_; ++y) {
          for (x = 0; x < tiff_width_; ++x) {
            if (*++p >= 0xff000000) {
              left = x;
              top = y;
              y = tiff_u_height_;
              break;
            }
          }
        }

        p = ((uint32_t*)data) +
            static_cast<ptrdiff_t>(tiff_width_) * tiff_u_height_;
        for (y = tiff_u_height_ - 1; y >= 0; --y) {
          for (x = tiff_width_ - 1; x >= 0; --x) {
            if (*--p >= 0xff000000) {
              right = x;
              bottom = y;
              y = -1;
              break;
            }
          }
        }

        if (right < left) {
          std::swap(left, right);
        }

        p = ((uint32_t*)data) + static_cast<ptrdiff_t>(tiff_width_) * top;
        for (y = top; y <= bottom; ++y) {
          for (x = 0; x < left; ++x) {
            if (p[x] >= 0xff000000) {
              left = x;
              break;
            }
          }

          for (x = tiff_width_ - 1; x > right; --x) {
            if (p[x] >= 0xff000000) {
              right = x;
              break;
            }
          }

          p += tiff_width_;
        }
      } break;
      case 16: {
        uint16_t* p = ((uint16_t*)data) + 3;

        for (y = 0; y < tiff_u_height_; ++y) {
          for (x = 0; x < tiff_width_; ++x) {
            if (p[x << 2] == 0xffff) {
              left = x;
              top = y;
              y = tiff_u_height_;
              break;
            }
          }

          p += tiff_width_ << 2;
        }

        p = ((uint16_t*)data) +
            (static_cast<ptrdiff_t>(tiff_width_ << 2)) * (tiff_u_height_ - 1) +
            3;
        for (y = tiff_u_height_ - 1; y >= 0; --y) {
          for (x = tiff_width_ - 1; x >= 0; --x) {
            if (p[x << 2] == 0xffff) {
              right = x;
              bottom = y;
              y = -1;
              break;
            }
          }

          p -= tiff_width_ << 2;
        }

        if (right < left) {
          std::swap(left, right);
        }

        p = ((uint16_t*)data) +
            (static_cast<ptrdiff_t>(tiff_width_ << 2)) * top + 3;
        for (y = top; y <= bottom; ++y) {
          for (x = 0; x < left; ++x) {
            if (p[x << 2] == 0xffff) {
              left = x;
              break;
            }
          }

          for (x = tiff_width_ - 1; x > right; --x) {
            if (p[x << 2] == 0xffff) {
              right = x;
              break;
            }
          }

          p += tiff_width_ << 2;
        }

      } break;
    }

    width_ = right + 1 - left;
    height_ = bottom + 1 - top;

    xpos_ += left;
    ypos_ += top;

    /***********************************************************************
     * Inpaint
     ***********************************************************************/
    int temp_copy;
    uint32_t a;
    uint32_t b;
    uint32_t c;
    uint32_t d;
    uint32_t* this_line = nullptr;
    uint32_t* prev_line = nullptr;
    auto* threadpool = mt::GetInstance();

    tiff_mask_ = std::make_unique<utils::Flex>(width_, height_);
    auto dt = utils::Flex(width_, height_);
    int mc;

    int n_threads = (std::max)(2, threadpool->GetNThreads());

    auto thread_lines = std::vector<std::vector<uint32_t>>(n_threads);
    auto thread_comp_lines = std::vector<std::vector<uint32_t>>(n_threads);

    for (int i = 0; i < n_threads; ++i) {
      thread_lines[i].resize(width_);
      thread_comp_lines[i].resize(width_);
    }

    uint32_t* bitmap32 = nullptr;
    uint64_t* bitmap64 = nullptr;
    if (bpp_ == 8) {
      bitmap32 =
          ((uint32_t*)data) + static_cast<ptrdiff_t>(top) * tiff_width_ + left;
    } else {
      bitmap64 =
          ((uint64_t*)data) + static_cast<ptrdiff_t>(top) * tiff_width_ + left;
    }

    auto tasks = mt::MultiFuture<std::pair<uint32_t*, int>>{};
    for (y = 0; y < height_; ++y) {
      int t = y % n_threads;
      this_line = thread_lines[t].data();
      uint32_t* comp = thread_comp_lines[t].data();
      bool first;

      x = 0;

      if (tasks.size() == n_threads) {
        tasks.wait();
        for (auto [comp_line, length] : tasks.get()) {
          dt.Copy((uint8_t*)comp_line, length);
          dt.NextLine();
        }
        tasks = {};
      }

      while (x < width_) {
        mc = 0;
        first = true;

        while (x < width_ && (bpp_ == 8 ? (bitmap32[x] < 0xff000000)
                                        : (bitmap64[x] < 0xffff000000000000))) {
          if (y == 0) {
            if (x == 0) {
              this_line[0] = 0x80000000;
            } else {
              this_line[x] = this_line[x - 1] + 3;
              if (bpp_ == 8) {
                bitmap32[x] = bitmap32[x - 1];
              } else {
                bitmap64[x] = bitmap64[x - 1];
              }
            }
          } else {
            if (x == 0) {
              b = prev_line[x] + 3;
              c = prev_line[x + 1] + 4;

              if (b < c) {
                d = b;
                temp_copy = x - tiff_width_;
              } else {
                d = c;
                temp_copy = x - tiff_width_ + 1;
              }
              d = b < c ? b : c;

              this_line[x] = d;
              if (bpp_ == 8) {
                bitmap32[x] = bitmap32[temp_copy];
              } else {
                bitmap64[x] = bitmap64[temp_copy];
              }

              a = b + 1;
              b = c - 1;
              d += 3;
            } else {
              if (first) {
                a = prev_line[x - 1] + 4;
                b = prev_line[x] + 3;
                d = this_line[x - 1] + 3;
                first = false;
              }

              if (x < width_ - 1) {
                c = prev_line[x + 1] + 4;
              } else {
                c = 0xffffffff;
              }
              temp_copy = x - 1;

              if (a < d) {
                d = a;
                temp_copy = x - tiff_width_ - 1;
              }
              if (b < d) {
                d = b;
                temp_copy = x - tiff_width_;
              }
              if (c < d) {
                d = c;
                temp_copy = x - tiff_width_ + 1;
              }

              this_line[x] = d;
              if (bpp_ == 8) {
                bitmap32[x] = bitmap32[temp_copy];
              } else {
                bitmap64[x] = bitmap64[temp_copy];
              }

              a = b + 1;
              b = c - 1;
              d += 3;
            }
          }
          ++x;
          ++mc;
        }

        if (mc != 0) {
          tiff_mask_->Write32(mc);
          mc = 0;
        }

        switch (bpp_) {
          case 8: {
            while (x < width_ && (bitmap32[x] >= 0xff000000)) {
              this_line[x++] = 0;
              ++mc;
            }
          } break;
          case 16: {
            while (x < width_ && (bitmap64[x] >= 0xffff000000000000)) {
              this_line[x++] = 0;
              ++mc;
            }
            //      }
          } break;
        }

        if (mc != 0) {
          tiff_mask_->Write32(0x80000000 | mc);
          mc = 0;
        }
      }

      if (bpp_ == 8) {
        bitmap32 += tiff_width_;
      } else {
        bitmap64 += tiff_width_;
      }

      if (y < height_ - 1) {
        tasks.push_back(threadpool->Queue([=, this] {
          int p = utils::CompressDTLine(this_line, (uint8_t*)comp, width_);
          return std::pair{comp, p};
        }));
      }

      tiff_mask_->NextLine();
      prev_line = this_line;
    }

    tasks.wait();
    for (auto [comp_line, length] : tasks.get()) {
      dt.Copy((uint8_t*)comp_line, length);
      dt.NextLine();
    }

    // backward
    int current_count = 0;
    int current_step;
    uint32_t dt_val;

    uint32_t mask;

    prev_line = thread_lines[(y - 2) % n_threads].data();

    // first line
    x = width_ - 1;
    if (bpp_ == 8) {
      bitmap32 -= tiff_width_;
    } else {
      bitmap64 -= tiff_width_;
    }

    while (x >= 0) {
      mask = tiff_mask_->ReadBackwards32();
      if ((mask & 0x80000000) != 0u) {  // solid
        x -= mask & 0x7fffffff;
        d = 3;
      } else {
        if (x == width_ - 1) {
          d = this_line[x] + 3;
          --mask;
          --x;
        }
        while (mask != 0u) {
          uint32_t best = this_line[x];
          if (d < best) {
            this_line[x] = d;
            if (bpp_ == 8) {
              bitmap32[x] = bitmap32[x + 1];
            } else {
              bitmap64[x] = bitmap64[x + 1];
            }
            d += 3;
          } else {
            d = best + 3;
          }
          --x;
          --mask;
        }
      }
    }

    std::swap(this_line, prev_line);

    // other lines
    for (y = height_ - 2; y >= 0; --y) {
      x = width_ - 1;
      c = d = 0x80000000;
      if (bpp_ == 8) {
        bitmap32 -= tiff_width_;
      } else {
        bitmap64 -= tiff_width_;
      }

      while (x >= 0) {
        mask = tiff_mask_->ReadBackwards32();

        if ((mask & 0x80000000) != 0u) {  // solid
          mc = mask & 0x7fffffff;
          x -= mc;
          memset(&this_line[x + 1], 0, mc << 2);
          d = 3;
        } else {
          b = prev_line[x] + 3;
          c = (x < width_ - 1) ? prev_line[x + 1] + 4 : 0x80000000;
          while (mask != 0u) {
            a = (x > 0) ? prev_line[x - 1] + 4 : 0x80000000;
            utils::ReadInpaintDT(dt, current_count, current_step, dt_val);
            int copy = 0;
            uint32_t best = dt_val;
            if (a < best) {
              best = a;
              copy = tiff_width_ + x - 1;
            }
            if (b < best) {
              best = b;
              copy = tiff_width_ + x;
            }
            if (c < best) {
              best = c;
              copy = tiff_width_ + x + 1;
            }
            if (d < best) {
              best = d;
              copy = x + 1;
            }

            if (copy != 0) {
              if (bpp_ == 8) {
                bitmap32[x] = bitmap32[copy];
              } else {
                bitmap64[x] = bitmap64[copy];
              }
            }

            this_line[x--] = best;
            c = b + 1;
            d = best + 3;
            b = a - 1;
            mask--;
          }
        }
      }

      std::swap(this_line, prev_line);
    }
  } else {
    width_ = tiff_width_;
    height_ = tiff_height_;

    tiff_mask_ = std::make_unique<utils::Flex>(width_, height_);

    for (y = 0; y < height_; ++y) {
      tiff_mask_->Write32(0x80000000 | width_);
      tiff_mask_->NextLine();
    }
  }

  /***********************************************************************
   * Extract channels
   ***********************************************************************/
  std::size_t channel_bytes = ((std::size_t)width_ * height_) << (bpp_ >> 4);

  for (int c = 0; c < 3; ++c) {
    channels_.emplace_back(channel_bytes);
  }

  if (spp_ == 4) {
    switch (bpp_) {
      case 8: {
        uint32_t* line = ((uint32_t*)data) +
                         static_cast<ptrdiff_t>(top) * tiff_width_ + left;
        int p = 0;
        for (y = 0; y < height_; ++y) {
          for (x = 0; x < width_; ++x) {
            uint32_t pixel = line[x];
            ((uint8_t*)channels_[0].data_.get())[p + x] = pixel & 0xff;
            ((uint8_t*)channels_[1].data_.get())[p + x] = (pixel >> 8) & 0xff;
            ((uint8_t*)channels_[2].data_.get())[p + x] = (pixel >> 16) & 0xff;
          }
          p += width_;
          line += tiff_width_;
        }
      } break;
      case 16: {
        uint64_t* line = ((uint64_t*)data) +
                         static_cast<ptrdiff_t>(top) * tiff_width_ + left;
        int p = 0;
        for (y = 0; y < height_; ++y) {
          for (x = 0; x < width_; ++x) {
            uint64_t pixel = line[x];
            ((uint16_t*)channels_[0].data_.get())[p + x] = pixel & 0xffff;
            ((uint16_t*)channels_[1].data_.get())[p + x] =
                (pixel >> 16) & 0xffff;
            ((uint16_t*)channels_[2].data_.get())[p + x] =
                (pixel >> 32) & 0xffff;
          }
          p += width_;
          line += tiff_width_;
        }
      } break;
    }
  } else {
    switch (bpp_) {
      case 8: {
        uint8_t byte;
        auto* bytes = (uint8_t*)data;
        std::size_t p = 0;
        for (y = 0; y < height_; ++y) {
          for (x = 0; x < width_; ++x) {
            byte = *bytes++;
            ((uint8_t*)channels_[0].data_.get())[p] = byte;
            // channel_totals[0] += gamma ? byte*byte : byte;

            byte = *bytes++;
            ((uint8_t*)channels_[1].data_.get())[p] = byte;
            // channel_totals[1] += gamma ? byte * byte : byte;

            byte = *bytes++;
            ((uint8_t*)channels_[2].data_.get())[p] = byte;
            // channel_totals[2] += gamma ? byte * byte : byte;

            ++p;
          }
        }
      } break;
      case 16: {
        uint16_t word;
        auto* words = (uint16_t*)data;
        std::size_t p = 0;
        for (y = 0; y < height_; ++y) {
          for (x = 0; x < width_; ++x) {
            word = *words++;
            ((uint16_t*)channels_[0].data_.get())[p] = word;
            // channel_totals[0] += gamma ? word * word : word;

            word = *words++;
            ((uint16_t*)channels_[1].data_.get())[p] = word;
            // channel_totals[1] += gamma ? word * word : word;

            word = *words++;
            ((uint16_t*)channels_[2].data_.get())[p] = word;
            // channel_totals[2] += gamma ? word * word : word;

            ++p;
          }
        }
      } break;
    }

    //  total_pixels += width * height;
  }
}

/***********************************************************************
 * Debugging
 ***********************************************************************/
void Image::MaskPng(int i) {
  int width = masks_[0].width_;
  int height = masks_[0].height_;  // +1 + masks[1]->height;

  std::size_t size = (std::size_t)width * height;
  auto temp = std::make_unique<uint8_t[]>(size);
  memset(temp.get(), 32, size);

  int px = 0;
  int py = 0;
  uint32_t cur;

  for (int l = 0; l < (int)masks_.size(); ++l) {
    auto* data = (uint32_t*)masks_[l].data_.get();
    uint8_t* line =
        temp.get() + static_cast<ptrdiff_t>(py) * masks_[0].width_ + px;
    for (int y = 0; y < masks_[l].height_; ++y) {
      int x = 0;
      while (x < masks_[l].width_) {
        float val;
        int count;

        cur = *data++;

        if ((cur & 0x80000000) != 0u) {
          count = cur & 0x00ffffff;
          if ((cur & 0x20000000) != 0u) {
            val = *(float*)data++;
          } else {
            val = (float)((cur >> 30) & 1);
          }
        } else {
          val = *((float*)&cur);
          count = 1;
        }

        int t = x + count;
        while (x < t) {
          line[x++] = (uint8_t)(val * 255);
        }
      }

      line += masks_[0].width_;
    }
    if ((l & 1) != 0) {
      px += masks_[l].width_ + 1;
    } else {
      py += masks_[l].height_ + 1;
    }
    break;
  }

  std::string filename = fmt::format("masks-{}.png", i);
  io::png::Pnger::Quick(filename.c_str(), temp.get(), width, height, width,
                        io::png::ColorType::GRAY);
}

}  // namespace multiblend::io
