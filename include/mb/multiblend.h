#pragma once

#include <array>

#include "mb/functions.h"
#include "mb/image.h"
#include "mb/mapalloc.h"

namespace multiblend {

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
  io::ImageType output_type = io::ImageType::MB_NONE;
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

  char* seamsave_filename = nullptr;
  char* seamload_filename = nullptr;
  char* xor_filename = nullptr;
};

struct Result {
  int output_bpp = 0;
  int width = 0;
  int height = 0;
  bool no_mask = false;

  std::array<memory::AlignedAllocPtr<void>, 3> output_channels;
  int min_xpos = 0x7fffffff;
  int min_ypos = 0x7fffffff;

  utils::Flex full_mask;
  TimingResult timing = {};
};

Result Multiblend(std::vector<io::Image>& images, Options opts);

}  // namespace multiblend
