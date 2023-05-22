#pragma once

#include <memory>

#include <tiffio.h>

namespace multiblend::io::tiff {
class CloseDeleter {
 public:
  void operator()(TIFF* tiff) const noexcept { TIFFClose(tiff); }
};

class FreeDeleter {
 public:
  void operator()(void* data) const noexcept { _TIFFfree(data); }
};

using TiffPtr = std::unique_ptr<TIFF, CloseDeleter>;

}  // namespace multiblend::io::tiff
