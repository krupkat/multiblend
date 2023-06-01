// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

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

struct GeoTIFFInfo {
  double XGeoRef, YGeoRef;
  double XCellRes, YCellRes;
  double projection[16];
  int nodata;
  bool set;
};

int geotiff_read(TIFF* tiff, GeoTIFFInfo* info);

int geotiff_write(TIFF* tiff, GeoTIFFInfo* info);

}  // namespace multiblend::io::tiff
