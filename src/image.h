#pragma once

#include <cstdio>

#include <jpeglib.h>
#include <png.h>
#include <tiffio.h>

#include "src/functions.h"
#include "src/geotiff.h"
#include "src/mapalloc.h"
#include "src/pyramid.h"

namespace multiblend::io {

enum class ImageType { MB_NONE, MB_TIFF, MB_JPEG, MB_PNG };

class Channel {
 public:
  explicit Channel(std::size_t bytes) : bytes_(bytes) {
    data_ = memory::MapAllocPtr{memory::MapAlloc::Alloc(bytes_),
                                memory::MapAllocDeleter{}};
  };

  memory::MapAllocPtr data_ = nullptr;
  std::size_t bytes_;
  FILE* file_ = nullptr;
};

class Image {
 public:
  explicit Image(char* filename);
  ~Image();

  Image(const Image&) = delete;
  Image& operator=(const Image&) = delete;
  Image(Image&&) = default;
  Image& operator=(Image&&) = default;

  char* filename_;
  ImageType type_;
  int width_;
  int height_;
  int xpos_;
  int ypos_;
  int xpos_add_ = 0;
  int ypos_add_ = 0;
  std::vector<Channel> channels_;
  Pyramid* pyramid_ = nullptr;
  tiff::GeoTIFFInfo geotiff_;
  int tiff_width_;
  int tiff_height_;
  int tiff_u_height_;
  int rows_per_strip_;
  int first_strip_;
  int end_strip_;
  uint16_t bpp_;
  uint16_t spp_;

  void Open();

  void Read(void* data, bool gamma);

  // std::size_t untrimmed_pixels;
  std::size_t untrimmed_bytes_;
  utils::Flex* tiff_mask_;
  float tiff_xres_, tiff_yres_;
  uint64_t mask_state_;
  int mask_count_;
  int mask_limit_;
  bool seam_present_;
  std::vector<utils::Flex*> masks_;

  void MaskPng(int i);

 private:
  TIFF* tiff_;
  FILE* file_;
  struct jpeg_decompress_struct cinfo_;
  struct jpeg_error_mgr jerr_;
  png_structp png_ptr_;
};

}  // namespace multiblend::io
