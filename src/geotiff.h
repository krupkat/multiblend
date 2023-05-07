#pragma once

#include <tiffio.h>

namespace multiblend::io::tiff {

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
