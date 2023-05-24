#pragma once

#include <cstdio>
#include <memory>
#include <optional>

#include "src/file.h"
#include "src/functions.h"
#ifdef MULTIBLEND_WITH_JPEG
#include "src/jpeg.h"
#endif
#include "src/mapalloc.h"
#include "src/pnger.h"
#include "src/pyramid.h"
#ifdef MULTIBLEND_WITH_TIFF
#include "src/tiff.h"
#endif

namespace multiblend::io {

enum class ImageType { MB_NONE, MB_TIFF, MB_JPEG, MB_PNG };

class Channel {
 public:
  explicit Channel(std::size_t bytes) : bytes_(bytes) {
    data_ = memory::MapAllocPtr<void>{memory::MapAlloc::Alloc(bytes_)};
  };

  memory::MapAllocPtr<void> data_ = nullptr;
  std::size_t bytes_;
};

class Image {
 public:
  explicit Image(char* filename);

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
  std::optional<Pyramid> pyramid_;
#ifdef MULTIBLEND_WITH_TIFF
  tiff::GeoTIFFInfo geotiff_;
#endif
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
  std::unique_ptr<utils::Flex> tiff_mask_;
  float tiff_xres_;
  float tiff_yres_;
  uint64_t mask_state_;
  int mask_count_;
  int mask_limit_;
  bool seam_present_;
  std::vector<utils::Flex> masks_;

  void MaskPng(int i);

 private:
  FilePtr file_;

#ifdef MULTIBLEND_WITH_TIFF
  tiff::TiffPtr tiff_;
#endif
#ifdef MULTIBLEND_WITH_JPEG
  std::unique_ptr<jpeg_error_mgr> jerr_;
  std::unique_ptr<jpeg_decompress_struct, jpeg::DecompressDeleter> cinfo_;
#endif
#ifdef MULTIBLEND_WITH_PNG
  std::unique_ptr<png_struct, png::PngReadStructDeleter> png_ptr_;
#endif
};

}  // namespace multiblend::io
