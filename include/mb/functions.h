// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <chrono>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include "mb/flex.h"
#include "mb/pyramid.h"

namespace multiblend::utils {

/***********************************************************************
 * Timer
 ***********************************************************************/
class Timer {
 public:
  void Start() { start_time_ = std::chrono::high_resolution_clock::now(); };

  double Read() {
    std::chrono::duration<double> elapsed =
        std::chrono::high_resolution_clock::now() - start_time_;
    return elapsed.count();
  };

  void Report(const char* name);

 private:
  std::chrono::high_resolution_clock::time_point start_time_;
};

void ShrinkMasks(std::vector<Flex>& masks, int n_levels);

void CompositeLine(const float* input_p, float* output_p, int i, int x_offset,
                   int in_level_width, int out_level_width, int out_level_pitch,
                   uint8_t* _mask, std::size_t mask_p);

void ReadInpaintDT(Flex& flex, int& current_count, int& current_step,
                   uint32_t& dt_val);

void ReadSeamDT(Flex& flex, int& current_count, int64_t& current_step,
                uint64_t& dt_val);

int CompressDTLine(const uint32_t* input, uint8_t* output, int width);

int CompressSeamLine(const uint64_t* input, uint8_t* output, int width);

void SwapH(Pyramid* py);

void UnswapH(Pyramid* py);

void SwapV(Pyramid* py);

void UnswapV(Pyramid* py);

}  // namespace multiblend::utils
