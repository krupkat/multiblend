#pragma once

#include "src/functions.h"
#include "src/image.h"

struct TimingResult {
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
};

struct Options {
  ImageType output_type = ImageType::MB_NONE;
  int output_bpp = 0;
  int fixed_levels = 0;
  bool wideblend = false;
  int add_levels = 0;
  bool all_threads = true;
  bool reverse = false;
  int wrap = 0;
  bool dither = true;
  bool gamma = false;
  bool no_mask = false;

  char* seamsave_filename = NULL;
  char* seamload_filename = NULL;
  char* xor_filename = NULL;
};

struct Result {
  int output_bpp = 0;
  int width = 0;
  int height = 0;
  bool no_mask = false;

  std::array<void*, 3> output_channels = {NULL, NULL, NULL};
  int min_xpos = 0x7fffffff;
  int min_ypos = 0x7fffffff;

  Flex full_mask;
  TimingResult timing = {};
};

Result DoWork(std::vector<Image*>& images, Options opts);