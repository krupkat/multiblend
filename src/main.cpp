/*
        Multiblend 2.0 (c) 2021 David Horman

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program. If not, see <https://www.gnu.org/licenses/>.

        The author can be contacted at davidhorman47@gmail.com
*/

#include <algorithm>
#include <array>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>

#include <jpeglib.h>
#include <tiffio.h>

extern int verbosity;

#include "src/functions.h"
#include "src/image.h"
#include "src/linux_overrides.h"
#include "src/mapalloc.h"
#include "src/multiblend.h"
#include "src/pnger.h"

int main(int argc, char* argv[]) {
  // This is here because of a weird problem encountered during development with
  // Visual Studio. It should never be triggered.
  if (verbosity != 1) {
    printf("bad compile?\n");
    exit(EXIT_FAILURE);
  }

  Timer timer_all, timer;
  timer_all.Start();

  TIFFSetWarningHandler(NULL);

  /***********************************************************************
   * Variables
   ***********************************************************************/
  std::vector<Image*> images;
  int fixed_levels = 0;
  int add_levels = 0;

  int width = 0;
  int height = 0;

  bool no_mask = false;
  bool big_tiff = false;
  bool bgr = false;
  bool wideblend = false;
  bool reverse = false;
  bool timing = false;
  bool dither = true;
  bool gamma = false;
  bool all_threads = true;
  int wrap = 0;

  TIFF* tiff_file = NULL;
  FILE* jpeg_file = NULL;
  Pnger* png_file = NULL;
  ImageType output_type = ImageType::MB_NONE;
  int jpeg_quality = -1;
  int compression = -1;
  char* seamsave_filename = NULL;
  char* seamload_filename = NULL;
  char* xor_filename = NULL;
  char* output_filename = NULL;
  int output_bpp = 0;

  double images_time = 0;
  double copy_time = 0;
  double seam_time = 0;
  double shrink_mask_time = 0;
  double shrink_time = 0;
  double laplace_time = 0;
  double blend_time = 0;
  double collapse_time = 0;
  double wrap_time = 0;
  double out_time = 0;
  double write_time = 0;

  /***********************************************************************
   * Help
   ***********************************************************************/
  if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help") ||
      !strcmp(argv[1], "/?")) {
    Output(1, "\n");
    Output(1,
           "Multiblend v2.0.0 (c) 2021 David Horman        "
           "http://horman.net/multiblend/\n");
    Output(1,
           "-------------------------------------------------------------------"
           "---------\n");

    printf(
        "Usage: multiblend [options] [-o OUTPUT] INPUT [X,Y] [INPUT] [X,Y] "
        "[INPUT]...\n");
    printf("\n");
    printf("Options:\n");
    printf("  --levels X / -l X      X: set number of blending levels to X\n");
    printf(
        "                        -X: decrease number of blending levels by "
        "X\n");
    printf(
        "                        +X: increase number of blending levels by "
        "X\n");
    printf(
        "  --depth D / -d D       Override automatic output image depth (8 or "
        "16)\n");
    printf("  --bgr                  Swap RGB order\n");
    printf(
        "  --wideblend            Calculate number of levels based on output "
        "image size,\n");
    printf("                         rather than input image size\n");
    printf(
        "  -w, --wrap=[mode]      Blend around images boundaries (NONE "
        "(default),\n");
    printf(
        "                         HORIZONTAL, VERTICAL). When specified "
        "without a mode,\n");
    printf("                         defaults to HORIZONTAL.\n");
    printf(
        "  --compression=X        Output file compression. For TIFF output, X "
        "may be:\n");
    printf("                         NONE (default), PACKBITS, or LZW\n");
    printf(
        "                         For JPEG output, X is JPEG quality (0-100, "
        "default 75)\n");
    printf(
        "                         For PNG output, X is PNG filter (0-9, "
        "default 3)\n");
    printf(
        "  --cache-threshold=     Allocate memory beyond X "
        "bytes/[K]ilobytes/\n");
    printf("      X[K/M/G]           [M]egabytes/[G]igabytes to disk\n");
    printf("  --no-dither            Disable dithering\n");
    printf(
        "  --tempdir <dir>        Specify temporary directory (default: system "
        "temp)\n");
    printf(
        "  --save-seams <file>    Save seams to PNG file for external "
        "editing\n");
    printf("  --load-seams <file>    Load seams from PNG file\n");
    printf(
        "  --no-output            Do not blend (for use with --save-seams)\n");
    printf(
        "                         Must be specified as last option before "
        "input images\n");
    printf("  --bigtiff              BigTIFF output\n");
    printf(
        "  --reverse              Reverse image priority (last=highest) for "
        "resolving\n");
    printf("                         indeterminate pixels\n");
    printf("  --quiet                Suppress output (except warnings)\n");
    printf("  --all-threads          Use all available CPU threads\n");
    printf(
        "  [X,Y]                  Optional position adjustment for previous "
        "input image\n");
    exit(EXIT_SUCCESS);
  }

  /***********************************************************************
  ************************************************************************
  * Parse arguments
  ************************************************************************
  ***********************************************************************/
  std::vector<char*> my_argv;

  bool skip = false;

  for (int i = 1; i < argc; ++i) {
    my_argv.push_back(argv[i]);

    if (!skip) {
      int c = 0;

      while (argv[i][c]) {
        if (argv[i][c] == '=') {
          argv[i][c++] = 0;
          if (argv[i][c]) {
            my_argv.push_back(&argv[i][c]);
          }
          break;
        }
        ++c;
      }

      if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--output")) {
        skip = true;
      }
    }
  }

  if ((int)my_argv.size() < 3)
    die("Error: Not enough arguments (try -h for help)");

  int pos;
  for (pos = 0; pos < (int)my_argv.size(); ++pos) {
    if (!strcmp(my_argv[pos], "-d") || !strcmp(my_argv[pos], "--d") ||
        !strcmp(my_argv[pos], "--depth") || !strcmp(my_argv[pos], "--bpp")) {
      if (++pos < (int)my_argv.size()) {
        output_bpp = atoi(my_argv[pos]);
        if (output_bpp != 8 && output_bpp != 16) {
          die("Error: Invalid output depth specified");
        }
      } else {
        die("Error: Missing parameter value");
      }
    } else if (!strcmp(my_argv[pos], "-l") ||
               !strcmp(my_argv[pos], "--levels")) {
      if (++pos < (int)my_argv.size()) {
        int n;
        if (my_argv[pos][0] == '+' || my_argv[pos][0] == '-') {
          sscanf_s(my_argv[pos], "%d%n", &add_levels, &n);
        } else {
          sscanf_s(my_argv[pos], "%d%n", &fixed_levels, &n);
          if (fixed_levels == 0) fixed_levels = 1;
        }
        if (my_argv[pos][n]) die("Error: Bad --levels parameter");
      } else {
        die("Error: Missing parameter value");
      }
    } else if (!strcmp(my_argv[pos], "--wrap") || !strcmp(my_argv[pos], "-w")) {
      if (pos + 1 >= (int)my_argv.size()) {
        die("Error: Missing parameters");
      }
      if (!strcmp(my_argv[pos + 1], "none") ||
          !strcmp(my_argv[pos + 1], "open"))
        ++pos;
      else if (!strcmp(my_argv[pos + 1], "horizontal") ||
               !strcmp(my_argv[pos + 1], "h")) {
        wrap = 1;
        ++pos;
      } else if (!strcmp(my_argv[pos + 1], "vertical") ||
                 !strcmp(my_argv[pos + 1], "v")) {
        wrap = 2;
        ++pos;
      } else if (!strcmp(my_argv[pos + 1], "both") ||
                 !strcmp(my_argv[pos + 1], "hv")) {
        wrap = 3;
        ++pos;
      } else
        wrap = 1;
    } else if (!strcmp(my_argv[pos], "--cache-threshold")) {
      if (pos + 1 >= (int)my_argv.size()) {
        die("Error: Missing parameters");
      }
      ++pos;
      int shift = 0;
      int n = 0;
      size_t len = strlen(my_argv[pos]);
      size_t threshold;
      sscanf_s(my_argv[pos], "%zu%n", &threshold, &n);
      if (n != len) {
        if (n == len - 1) {
          switch (my_argv[pos][len - 1]) {
            case 'k':
            case 'K':
              shift = 10;
              break;
            case 'm':
            case 'M':
              shift = 20;
              break;
            case 'g':
            case 'G':
              shift = 30;
              break;
            default:
              die("Error: Bad --cache-threshold parameter");
          }
          threshold <<= shift;
        } else {
          die("Error: Bad --cache-threshold parameter");
        }
      }
      MapAlloc::CacheThreshold(threshold);
    } else if (!strcmp(my_argv[pos], "--nomask") ||
               !strcmp(my_argv[pos], "--no-mask"))
      no_mask = true;
    else if (!strcmp(my_argv[pos], "--timing") ||
             !strcmp(my_argv[pos], "--timings"))
      timing = true;
    else if (!strcmp(my_argv[pos], "--bigtiff"))
      big_tiff = true;
    else if (!strcmp(my_argv[pos], "--bgr"))
      bgr = true;
    else if (!strcmp(my_argv[pos], "--wideblend"))
      wideblend = true;
    else if (!strcmp(my_argv[pos], "--reverse"))
      reverse = true;
    else if (!strcmp(my_argv[pos], "--gamma"))
      gamma = true;
    else if (!strcmp(my_argv[pos], "--no-dither") ||
             !strcmp(my_argv[pos], "--nodither"))
      dither = false;
    //  else if (!strcmp(my_argv[i], "--force"))     force_coverage =
    // true;
    else if (!strncmp(my_argv[pos], "-f", 2))
      Output(0, "ignoring Enblend option -f\n");
    else if (!strcmp(my_argv[pos], "-a"))
      Output(0, "ignoring Enblend option -a\n");
    else if (!strcmp(my_argv[pos], "--no-ciecam"))
      Output(0, "ignoring Enblend option --no-ciecam\n");
    else if (!strcmp(my_argv[pos], "--primary-seam-generator")) {
      Output(0, "ignoring Enblend option --primary-seam-generator\n");
      ++pos;
    }

    else if (!strcmp(my_argv[pos], "--compression")) {
      if (++pos < (int)my_argv.size()) {
        if (strcmp(my_argv[pos], "0") == 0)
          jpeg_quality = 0;
        else if (atoi(my_argv[pos]) > 0)
          jpeg_quality = atoi(my_argv[pos]);
        else if (_stricmp(my_argv[pos], "lzw") == 0)
          compression = COMPRESSION_LZW;
        else if (_stricmp(my_argv[pos], "packbits") == 0)
          compression = COMPRESSION_PACKBITS;
        //    else if (_stricmp(my_argv[i], "deflate")
        //== 0) compression = COMPRESSION_DEFLATE;
        else if (_stricmp(my_argv[pos], "none") == 0)
          compression = COMPRESSION_NONE;
        else
          die("Error: Unknown compression codec %s", my_argv[pos]);
      } else {
        die("Error: Missing parameter value");
      }
    } else if (!strcmp(my_argv[pos], "-v") ||
               !strcmp(my_argv[pos], "--verbose"))
      ++verbosity;
    else if (!strcmp(my_argv[pos], "-q") || !strcmp(my_argv[pos], "--quiet"))
      --verbosity;
    else if ((!strcmp(my_argv[pos], "--saveseams") ||
              !strcmp(my_argv[pos], "--save-seams")) &&
             pos < (int)my_argv.size() - 1)
      seamsave_filename = my_argv[++pos];
    else if ((!strcmp(my_argv[pos], "--loadseams") ||
              !strcmp(my_argv[pos], "--load-seams")) &&
             pos < (int)my_argv.size() - 1)
      seamload_filename = my_argv[++pos];
    else if ((!strcmp(my_argv[pos], "--savexor") ||
              !strcmp(my_argv[pos], "--save-xor")) &&
             pos < (int)my_argv.size() - 1)
      xor_filename = my_argv[++pos];
    else if (!strcmp(my_argv[pos], "--tempdir") ||
             !strcmp(my_argv[pos], "--tmpdir") && pos < (int)my_argv.size() - 1)
      MapAlloc::SetTmpdir(my_argv[++pos]);
    else if (!strcmp(my_argv[pos], "--all-threads"))
      all_threads = true;
    else if (!strcmp(my_argv[pos], "-o") || !strcmp(my_argv[pos], "--output")) {
      if (++pos < (int)my_argv.size()) {
        output_filename = my_argv[pos];
        char* ext = strrchr(output_filename, '.');

        if (!ext) {
          die("Error: Unknown output filetype");
        }

        ++ext;
        if (!(_stricmp(ext, "jpg") && _stricmp(ext, "jpeg"))) {
          output_type = ImageType::MB_JPEG;
          if (jpeg_quality == -1) jpeg_quality = 75;
        } else if (!(_stricmp(ext, "tif") && _stricmp(ext, "tiff"))) {
          output_type = ImageType::MB_TIFF;
        } else if (!_stricmp(ext, "png")) {
          output_type = ImageType::MB_PNG;
        } else {
          die("Error: Unknown file extension");
        }

        ++pos;
        break;
      }
    } else if (!strcmp(my_argv[pos], "--no-output")) {
      ++pos;
      break;
    } else {
      die("Error: Unknown argument \"%s\"", my_argv[pos]);
    }
  }

  if (compression != -1) {
    if (output_type != ImageType::MB_TIFF) {
      Output(0,
             "Warning: non-TIFF output; ignoring TIFF compression setting\n");
    }
  } else if (output_type == ImageType::MB_TIFF) {
    compression = COMPRESSION_LZW;
  }

  if (jpeg_quality != -1 && output_type != ImageType::MB_JPEG &&
      output_type != ImageType::MB_PNG) {
    Output(0,
           "Warning: non-JPEG/PNG output; ignoring compression quality "
           "setting\n");
  }

  if ((jpeg_quality < -1 || jpeg_quality > 9) &&
      output_type == ImageType::MB_PNG) {
    die("Error: Bad PNG compression quality setting\n");
  }

  if (output_type == ImageType::MB_NONE && !seamsave_filename)
    die("Error: No output file specified");
  if (seamload_filename && seamsave_filename)
    die("Error: Cannot load and save seams at the same time");
  if (wrap == 3)
    die("Error: Wrapping in both directions is not currently supported");

  if (!strcmp(my_argv[pos], "--")) {
    ++pos;
  }

  /***********************************************************************
   * Push remaining arguments to images vector
   ***********************************************************************/

  while (pos < (int)my_argv.size()) {
    if (images.size()) {
      int x, y;
      int n = 0;
      sscanf_s(my_argv[pos], "%d,%d%n", &x, &y, &n);
      if (!my_argv[pos][n]) {
        images.back()->xpos_add = x;
        images.back()->ypos_add = y;
        pos++;
        continue;
      }
    }
    images.push_back(new Image(my_argv[pos++]));
  }

  int n_images = (int)images.size();

  if (n_images == 0) die("Error: No input files specified");
  if (seamsave_filename && n_images > 256) {
    seamsave_filename = NULL;
    Output(0, "Warning: seam saving not possible with more than 256 images");
  }
  if (seamload_filename && n_images > 256) {
    seamload_filename = NULL;
    Output(0, "Warning: seam loading not possible with more than 256 images");
  }
  if (xor_filename && n_images > 255) {
    xor_filename = NULL;
    Output(0, "Warning: XOR map saving not possible with more than 255 images");
  }

  /***********************************************************************
   * Print banner
   ***********************************************************************/
  Output(1, "\n");
  Output(1,
         "Multiblend v2.0.0 (c) 2021 David Horman        "
         "http://horman.net/multiblend/\n");
  Output(1,
         "---------------------------------------------------------------------"
         "-------\n");

  /***********************************************************************
  ************************************************************************
  * Open output
  ************************************************************************
  ***********************************************************************/
  switch (output_type) {
    case ImageType::MB_TIFF: {
      if (!big_tiff)
        tiff_file = TIFFOpen(output_filename, "w");
      else
        tiff_file = TIFFOpen(output_filename, "w8");
      if (!tiff_file) die("Error: Could not open output file");
    } break;
    case ImageType::MB_JPEG: {
      if (output_bpp == 16)
        die("Error: 16bpp output is incompatible with JPEG output");
      fopen_s(&jpeg_file, output_filename, "wb");
      if (!jpeg_file) die("Error: Could not open output file");
    } break;
    case ImageType::MB_PNG: {
      fopen_s(&jpeg_file, output_filename, "wb");
      if (!jpeg_file) die("Error: Could not open output file");
    } break;
  }

  /***********************************************************************
  ************************************************************************
  * Process images
  ************************************************************************
  ***********************************************************************/
  auto result = DoWork(images, {
                                   .output_type = output_type,
                                   .output_bpp = output_bpp,
                                   .fixed_levels = fixed_levels,
                                   .wideblend = wideblend,
                                   .add_levels = add_levels,
                                   .all_threads = all_threads,
                                   .reverse = reverse,
                                   .wrap = wrap,
                                   .dither = dither,
                                   .gamma = gamma,
                                   .no_mask = no_mask,
                                   .seamsave_filename = seamsave_filename,
                                   .seamload_filename = seamload_filename,
                                   .xor_filename = xor_filename,
                               });

/***********************************************************************
 * Write
 ***********************************************************************/
#define ROWS_PER_STRIP 64

  if (output_type != ImageType::MB_NONE) {
    Output(1, "Writing %s...\n", output_filename);

    timer.Start();

    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    JSAMPARRAY scanlines = NULL;

    int spp = result.no_mask ? 3 : 4;

    int bytes_per_pixel = spp << (result.output_bpp >> 4);
    int bytes_per_row = bytes_per_pixel * result.width;

    int n_strips = (int)((result.height + ROWS_PER_STRIP - 1) / ROWS_PER_STRIP);
    int remaining = result.height;
    void* strip =
        malloc((ROWS_PER_STRIP * (int64_t)result.width) * bytes_per_pixel);
    void* oc_p[3] = {result.output_channels[0], result.output_channels[1],
                     result.output_channels[2]};
    if (bgr) std::swap(oc_p[0], oc_p[2]);

    switch (output_type) {
      case ImageType::MB_TIFF: {
        TIFFSetField(tiff_file, TIFFTAG_IMAGEWIDTH, result.width);
        TIFFSetField(tiff_file, TIFFTAG_IMAGELENGTH, result.height);
        TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, compression);
        TIFFSetField(tiff_file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tiff_file, TIFFTAG_ROWSPERSTRIP, ROWS_PER_STRIP);
        TIFFSetField(tiff_file, TIFFTAG_BITSPERSAMPLE, result.output_bpp);
        if (result.no_mask) {
          TIFFSetField(tiff_file, TIFFTAG_SAMPLESPERPIXEL, 3);
        } else {
          TIFFSetField(tiff_file, TIFFTAG_SAMPLESPERPIXEL, 4);
          uint16_t out[1] = {EXTRASAMPLE_UNASSALPHA};
          TIFFSetField(tiff_file, TIFFTAG_EXTRASAMPLES, 1, &out);
        }

        TIFFSetField(tiff_file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        if (images[0]->tiff_xres != -1) {
          TIFFSetField(tiff_file, TIFFTAG_XRESOLUTION, images[0]->tiff_xres);
          TIFFSetField(tiff_file, TIFFTAG_XPOSITION,
                       (float)(result.min_xpos / images[0]->tiff_xres));
        }
        if (images[0]->tiff_yres != -1) {
          TIFFSetField(tiff_file, TIFFTAG_YRESOLUTION, images[0]->tiff_yres);
          TIFFSetField(tiff_file, TIFFTAG_YPOSITION,
                       (float)(result.min_ypos / images[0]->tiff_yres));
        }

        if (images[0]->geotiff.set) {
          // if we got a georeferenced input, store the geotags in the output
          GeoTIFFInfo info(images[0]->geotiff);
          info.XGeoRef = result.min_xpos * images[0]->geotiff.XCellRes;
          info.YGeoRef = -result.min_ypos * images[0]->geotiff.YCellRes;
          Output(1, "Output georef: UL: %f %f, pixel size: %f %f\n",
                 info.XGeoRef, info.YGeoRef, info.XCellRes, info.YCellRes);
          geotiff_write(tiff_file, &info);
        }
      } break;
      case ImageType::MB_JPEG: {
        cinfo.err = jpeg_std_error(&jerr);
        jpeg_create_compress(&cinfo);
        jpeg_stdio_dest(&cinfo, jpeg_file);

        cinfo.image_width = result.width;
        cinfo.image_height = result.height;
        cinfo.input_components = 3;
        cinfo.in_color_space = JCS_RGB;

        jpeg_set_defaults(&cinfo);
        jpeg_set_quality(&cinfo, jpeg_quality, true);
        jpeg_start_compress(&cinfo, true);
      } break;
      case ImageType::MB_PNG: {
        png_file = new Pnger(
            output_filename, NULL, result.width, result.height,
            result.no_mask ? PNG_COLOR_TYPE_RGB : PNG_COLOR_TYPE_RGB_ALPHA,
            result.output_bpp, jpeg_file, jpeg_quality);
      } break;
    }

    if (output_type == ImageType::MB_PNG || output_type == ImageType::MB_JPEG) {
      scanlines = new JSAMPROW[ROWS_PER_STRIP];
      for (int i = 0; i < ROWS_PER_STRIP; ++i) {
        scanlines[i] = (JSAMPROW) & ((uint8_t*)strip)[i * bytes_per_row];
      }
    }

    result.full_mask.Start();

    for (int s = 0; s < n_strips; ++s) {
      int strip_p = 0;
      int rows = std::min(remaining, ROWS_PER_STRIP);

      for (int strip_y = 0; strip_y < rows; ++strip_y) {
        int x = 0;
        while (x < result.width) {
          uint32_t cur = result.full_mask.ReadForwards32();
          if (cur & 0x80000000) {
            int lim = x + (cur & 0x7fffffff);
            switch (result.output_bpp) {
              case 8: {
                while (x < lim) {
                  ((uint8_t*)strip)[strip_p++] = ((uint8_t*)(oc_p[0]))[x];
                  ((uint8_t*)strip)[strip_p++] = ((uint8_t*)(oc_p[1]))[x];
                  ((uint8_t*)strip)[strip_p++] = ((uint8_t*)(oc_p[2]))[x];
                  if (!result.no_mask) ((uint8_t*)strip)[strip_p++] = 0xff;
                  ++x;
                }
              } break;
              case 16: {
                while (x < lim) {
                  ((uint16_t*)strip)[strip_p++] = ((uint16_t*)(oc_p[0]))[x];
                  ((uint16_t*)strip)[strip_p++] = ((uint16_t*)(oc_p[1]))[x];
                  ((uint16_t*)strip)[strip_p++] = ((uint16_t*)(oc_p[2]))[x];
                  if (!result.no_mask) ((uint16_t*)strip)[strip_p++] = 0xffff;
                  ++x;
                }
              } break;
            }
          } else {
            size_t t = (size_t)cur * bytes_per_pixel;
            switch (result.output_bpp) {
              case 8: {
                ZeroMemory(&((uint8_t*)strip)[strip_p], t);
              } break;
              case 16: {
                ZeroMemory(&((uint16_t*)strip)[strip_p], t);
              } break;
            }
            strip_p += cur * spp;
            x += cur;
          }
        }

        switch (result.output_bpp) {
          case 8: {
            oc_p[0] = &((uint8_t*)(oc_p[0]))[result.width];
            oc_p[1] = &((uint8_t*)(oc_p[1]))[result.width];
            oc_p[2] = &((uint8_t*)(oc_p[2]))[result.width];
          } break;
          case 16: {
            oc_p[0] = &((uint16_t*)(oc_p[0]))[result.width];
            oc_p[1] = &((uint16_t*)(oc_p[1]))[result.width];
            oc_p[2] = &((uint16_t*)(oc_p[2]))[result.width];
          } break;
        }
      }

      switch (output_type) {
        case ImageType::MB_TIFF: {
          TIFFWriteEncodedStrip(tiff_file, s, strip,
                                rows * (int64_t)bytes_per_row);
        } break;
        case ImageType::MB_JPEG: {
          jpeg_write_scanlines(&cinfo, scanlines, rows);
        } break;
        case ImageType::MB_PNG: {
          png_file->WriteRows(scanlines, rows);
        } break;
      }

      remaining -= ROWS_PER_STRIP;
    }

    switch (output_type) {
      case ImageType::MB_TIFF: {
        TIFFClose(tiff_file);
      } break;
      case ImageType::MB_JPEG: {
        jpeg_finish_compress(&cinfo);
        jpeg_destroy_compress(&cinfo);
        fclose(jpeg_file);
      } break;
    }

    write_time = timer.Read();
  }

  /***********************************************************************
   * Timing
   ***********************************************************************/
  if (timing) {
    printf("\n");
    printf("Images:   %.3fs\n", images_time);
    printf("Seaming:  %.3fs\n", seam_time);
    if (output_type != ImageType::MB_NONE) {
      printf("Masks:    %.3fs\n", shrink_mask_time);
      printf("Copy:     %.3fs\n", copy_time);
      printf("Shrink:   %.3fs\n", shrink_time);
      printf("Laplace:  %.3fs\n", laplace_time);
      printf("Blend:    %.3fs\n", blend_time);
      printf("Collapse: %.3fs\n", collapse_time);
      if (wrap) printf("Wrapping: %.3fs\n", wrap_time);
      printf("Output:   %.3fs\n", out_time);
      printf("Write:    %.3fs\n", write_time);
    }
  }

  /***********************************************************************
   * Clean up
   ***********************************************************************/
  if (timing) {
    if (output_type == ImageType::MB_NONE) {
      timer_all.Report("\nExecution complete. Total execution time");
    } else {
      timer_all.Report("\nBlend complete. Total execution time");
    }
  }

  return EXIT_SUCCESS;
}