#include "src/tiff.h"

#include <cstdint>
#include <cstdlib>

#include <tiffio.h>

namespace multiblend::io::tiff {

// some defintions for geotiff

constexpr int kTifftagGeoPixelScale = 33550;
constexpr int kTifftagGeoTiePoints = 33922;
constexpr int kTifftagGeoTransMatrix = 34264;
constexpr int kTifftagGeoKeyDirectory = 34735;
constexpr int kTifftagGeoDoubleParams = 34736;
constexpr int kTifftagGeoAsciiParams = 34737;

constexpr int kTifftagGdalNoData = 42113;

static const TIFFFieldInfo xtiffFieldInfo[] = {
    /* XXX Insert Your tags here */
    {kTifftagGeoPixelScale, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, 1, 1,
     (char*)"GeoPixelScale"},
    {kTifftagGeoTransMatrix, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, 1, 1,
     (char*)"GeoTransformationMatrix"},
    {kTifftagGeoTiePoints, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, 1, 1,
     (char*)"GeoTiePoints"},
    {kTifftagGeoKeyDirectory, -1, -1, TIFF_SHORT, FIELD_CUSTOM, 1, 1,
     (char*)"GeoKeyDirectory"},
    {kTifftagGeoDoubleParams, -1, -1, TIFF_DOUBLE, FIELD_CUSTOM, 1, 1,
     (char*)"GeoDoubleParams"},
    {kTifftagGeoAsciiParams, -1, -1, TIFF_ASCII, FIELD_CUSTOM, 1, 0,
     (char*)"GeoASCIIParams"},
    {kTifftagGdalNoData, -1, -1, TIFF_ASCII, FIELD_CUSTOM, 1, 0,
     (char*)"GDALNoDataValue"}};

void geotiff_register(TIFF* tif) {
  TIFFMergeFieldInfo(tif, xtiffFieldInfo,
                     sizeof(xtiffFieldInfo) / sizeof(xtiffFieldInfo[0]));
}

/** Read geotiff tags from an image. Only accept images with correct geocoding.

    Returns  1 if reading was successfull, 0 if it failed.
*/
int geotiff_read(TIFF* tiff, GeoTIFFInfo* info) {
  uint16_t nCount = 0;
  double* geo_scale;
  // clear geotiff info
  //  memset(info,0,sizeof(GeoTIFFInfo));
  geotiff_register(tiff);

  if ((TIFFGetField(tiff, kTifftagGeoPixelScale, &nCount, &geo_scale) == 0) ||
      nCount < 2) {
    return 0;
  }

  info->XCellRes = geo_scale[0];
  info->YCellRes = geo_scale[1];
  double* tiepoints;

  if ((TIFFGetField(tiff, kTifftagGeoTiePoints, &nCount, &tiepoints) == 0) ||
      nCount < 6) {
    return 0;
  }
  info->XGeoRef = tiepoints[3] - tiepoints[0] * (geo_scale[0]);
  info->YGeoRef = tiepoints[4] - tiepoints[1] * (geo_scale[1]);
  // TODO(dh): check if tiepoints refer to center of upper left pixel or upper
  // left edge of upper left pixel
  char* nodata;

  if (TIFFGetField(tiff, kTifftagGdalNoData, &nCount, &nodata) != 0) {
    info->nodata = std::atoi(nodata);
  }

  // TODO(dh): read coordinate system definitions...
  info->set = true;
  return 1;
}

/** Write geotiff tags to an image */
int geotiff_write(TIFF* tiff, GeoTIFFInfo* info) {
  geotiff_register(tiff);
  double scale[3];
  scale[0] = info->XCellRes;
  scale[1] = info->YCellRes;
  scale[2] = 0.0;
  TIFFSetField(tiff, kTifftagGeoPixelScale, 3, scale);
  double tiepoint[6];
  tiepoint[0] = tiepoint[1] = tiepoint[2] = tiepoint[5] = 0;
  tiepoint[3] = info->XGeoRef;
  tiepoint[4] = info->YGeoRef;
  TIFFSetField(tiff, kTifftagGeoTiePoints, 6, tiepoint);
  char nodata[50];
  //  SNPRINTF(nodata,50,"%d",info->nodata);
  TIFFSetField(tiff, kTifftagGdalNoData, nodata);
  return 1;
}

}  // namespace multiblend::io::tiff
