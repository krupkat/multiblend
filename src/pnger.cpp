#include "src/pnger.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include <png.h>

#include "src/functions.h"
#include "src/linux_overrides.h"

namespace multiblend::io::png {

std::vector<png_color> Pnger::palette_ = {};

Pnger::Pnger(const char* filename, const char* name, int width, int height,
             int type, int bpp, FilePtr file, int compression) {
  y_ = 0;
  height_ = height;

  if (type == PNG_COLOR_TYPE_PALETTE && palette_.empty()) {
    palette_.resize(256);

    double base = 2;
    double rad;
    double r;
    double g;
    double b;

    for (int i = 0; i < 255; ++i) {
      rad = base;
      r = std::max(0.0, std::min(1.0, std::min(rad, 4 - rad)));
      rad += 2;
      if (rad >= 6) {
        rad -= 6;
      }
      g = std::max(0.0, std::min(1.0, std::min(rad, 4 - rad)));
      rad += 2;
      if (rad >= 6) {
        rad -= 6;
      }
      b = std::max(0.0, std::min(1.0, std::min(rad, 4 - rad)));
      base += 6 * 0.618033988749895;
      if (base >= 6) {
        base -= 6;
      }
      palette_[i].red = (png_byte)std::lround(sqrt(r) * 255);
      palette_[i].green = (png_byte)std::lround(sqrt(g) * 255);
      palette_[i].blue = (png_byte)std::lround(sqrt(b) * 255);
    }

    palette_[255].red = 0;
    palette_[255].green = 0;
    palette_[255].blue = 0;
  }

  if (file == nullptr) {
    FILE* tmp_file = nullptr;
    fopen_s(&tmp_file, filename, "wb");
    if (tmp_file == nullptr) {
      if (name != nullptr) {
        utils::Output(0, "WARNING: Could not save %s\n", name);
      }
      return;
    }
    file_ = {tmp_file, FileDeleter{}};
  } else {
    file_ = std::move(file);
  }

  png_ptr_ = {
      png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr),
      PngWriteStructDeleter{}};
  if (png_ptr_ == nullptr) {
    if (name != nullptr) {
      utils::Output(0, "WARNING: PNG create failed\n");
    }
    file_.reset();
    remove(filename);
    return;
  }

  info_ptr_ = {png_create_info_struct(png_ptr_.get()),
               PngInfoStructDeleter{png_ptr_.get()}};
  if (info_ptr_ == nullptr) {
    if (name != nullptr) {
      utils::Output(0, "WARNING: PNG create failed\n");
    }
    png_ptr_.reset();
    file_.reset();
    remove(filename);
    return;
  }

  png_init_io(png_ptr_.get(), file_.get());

  png_set_IHDR(png_ptr_.get(), info_ptr_.get(), width, height, bpp, type,
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);
  if (type == PNG_COLOR_TYPE_PALETTE) {
    png_set_PLTE(png_ptr_.get(), info_ptr_.get(), palette_.data(), 256);
  }

  png_write_info(png_ptr_.get(), info_ptr_.get());
  png_set_compression_level(png_ptr_.get(), compression < 0 ? 3 : compression);
  if (bpp == 16) {
    png_set_swap(png_ptr_.get());
  }

  if (name != nullptr) {
    auto size = (type == PNG_COLOR_TYPE_RGB_ALPHA ? (width << 2)
                 : type == PNG_COLOR_TYPE_RGB     ? width * 3
                                                  : width)
                << (bpp >> 4);
    line_.resize(size);
  }
};

void Pnger::Write() {
  if (file_ == nullptr) {
    return;
  }

  png_write_row(png_ptr_.get(), line_.data());

  if (++y_ == height_) {
    //  printf("png close\n");
    png_write_end(png_ptr_.get(), nullptr);
    info_ptr_.reset();
    png_ptr_.reset();
    file_.reset();
  }
}

void Pnger::WriteRows(uint8_t** rows, int num_rows) {
  if (file_ == nullptr) {
    return;
  }

  png_write_rows(png_ptr_.get(), rows, num_rows);

  if ((y_ += num_rows) == height_) {
    png_write_end(png_ptr_.get(), nullptr);
    info_ptr_.reset();
    png_ptr_.reset();
    file_.reset();
  }
}

void Pnger::Quick(char* filename, uint8_t* data, int width, int height,
                  int pitch, int type) {
  Pnger temp(filename, nullptr, width, height, type);

  for (int y = 0; y < height; ++y) {
    temp.WriteRows(&data, 1);
    data += (type == PNG_COLOR_TYPE_RGB_ALPHA) ? (pitch << 2) : pitch;
  }
}

}  // namespace multiblend::io::png
