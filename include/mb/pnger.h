#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#ifdef MULTIBLEND_WITH_PNG
#include <png.h>
#endif

#include "mb/file.h"

namespace multiblend::io::png {

constexpr bool Enabled() {
#ifdef MULTIBLEND_WITH_PNG
  return true;
#else
  return false;
#endif
}

#ifdef MULTIBLEND_WITH_PNG
class PngWriteStructDeleter {
 public:
  void operator()(png_struct* png_ptr) const noexcept {
    png_destroy_write_struct(&png_ptr, nullptr);
  }
};

class PngReadStructDeleter {
 public:
  void operator()(png_struct* png_ptr) const noexcept {
    png_destroy_read_struct(&png_ptr, nullptr, nullptr);
  }
};

class PngInfoStructDeleter {
 public:
  PngInfoStructDeleter() = default;
  explicit PngInfoStructDeleter(png_struct* png_ptr) : png_ptr_(png_ptr) {}
  void operator()(png_info* info_ptr) const noexcept {
    if (png_ptr_ != nullptr) {
      png_destroy_info_struct(png_ptr_, &info_ptr);
    }
  }

  png_struct* png_ptr_ = nullptr;
};
#endif

enum class ColorType { RGB, RGB_ALPHA, PALETTE, GRAY };

class Pnger {
 public:
  Pnger(const char* filename, const char* name, int width, int height,
        ColorType type, int bpp = 8, FilePtr file = nullptr,
        int compression = -1);

  Pnger(const Pnger&) = delete;
  Pnger& operator=(const Pnger&) = delete;
  Pnger(Pnger&&) = default;
  Pnger& operator=(Pnger&&) = default;

  bool Ready() { return !(file_ == nullptr); };
  void WriteRows(uint8_t** rows, int num_rows);
  void Write();
  static void Quick(const char* filename, uint8_t* data, int width, int height,
                    int pitch, ColorType type);

  std::vector<uint8_t> line_;

 private:
  int y_;
  int height_;

  FilePtr file_;

#ifdef MULTIBLEND_WITH_PNG
  std::unique_ptr<png_struct, PngWriteStructDeleter> png_ptr_;
  std::unique_ptr<png_info, PngInfoStructDeleter> info_ptr_;

  static std::vector<png_color> palette_;
#endif
};

}  // namespace multiblend::io::png
