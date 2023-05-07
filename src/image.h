#pragma once

#include <stdio.h>

#include <jpeglib.h>
#include <png.h>
#include <tiffio.h>

#include "src/functions.h"
#include "src/geotiff.h"
#include "src/mapalloc.h"
#include "src/pyramid.h"

namespace multiblend {

enum class ImageType { MB_NONE, MB_TIFF, MB_JPEG, MB_PNG };

class Channel {
 public:
  Channel(size_t _bytes) : bytes(_bytes) { data = MapAlloc::Alloc(bytes); };
  ~Channel() { MapAlloc::Free(data); };
  void* data;
  size_t bytes;
  FILE* file = NULL;
};

class Image {
 public:
  Image(char* _filename);
  ~Image();
  char* filename;
  ImageType type;
  int width;
  int height;
  int xpos;
  int ypos;
  int xpos_add = 0;
  int ypos_add = 0;
  std::vector<Channel*> channels;
  Pyramid* pyramid = NULL;
  GeoTIFFInfo geotiff;
  int tiff_width;
  int tiff_height;
  int tiff_u_height;
  int rows_per_strip;
  int first_strip;
  int end_strip;
  uint16_t bpp;
  uint16_t spp;
  void Open();
  void Read(void* data, bool gamma);
  // size_t untrimmed_pixels;
  size_t untrimmed_bytes;
  utils::Flex* tiff_mask;
  float tiff_xres, tiff_yres;
  uint64_t mask_state;
  int mask_count;
  int mask_limit;
  bool seam_present;
  std::vector<utils::Flex*> masks;
  void MaskPng(int i);

 private:
  TIFF* tiff;
  FILE* file;
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  png_structp png_ptr;
};

}  // namespace multiblend
