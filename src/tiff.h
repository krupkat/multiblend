#pragma once

#include <cstdio>
#include <memory>

#include <tiffio.h>

namespace multiblend::io::tiff {
class TiffDeleter {
 public:
  void operator()(TIFF* tiff) const { TIFFClose(tiff); }
};

using TiffPtr = std::unique_ptr<TIFF, TiffDeleter>;

}  // namespace multiblend::io::tiff
