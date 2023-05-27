#pragma once

#include <memory>

#include <jpeglib.h>

namespace multiblend::io::jpeg {

class DecompressDeleter {
 public:
  void operator()(jpeg_decompress_struct* cinfo) const noexcept {
    jpeg_finish_decompress(cinfo);
    jpeg_destroy_decompress(cinfo);
    delete cinfo;
  }
};

class CompressDeleter {
 public:
  void operator()(jpeg_compress_struct* cinfo) const noexcept {
    jpeg_finish_compress(cinfo);
    jpeg_destroy_compress(cinfo);
    delete cinfo;
  }
};

}  // namespace multiblend::io::jpeg
