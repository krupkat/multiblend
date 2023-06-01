// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <chrono>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include "mb/pyramid.h"

namespace multiblend::utils {

class FreeDeleter {
 public:
  void operator()(void* ptr) const { free(ptr); }
};

/***********************************************************************
 * Flexible data class
 ***********************************************************************/
class Flex {
 public:
  Flex() : width_(0), height_(0) {}
  Flex(int width, int height) : width_(width), height_(height) {
    rows_.resize(height_);
    NextLine();
  };

  Flex(const Flex& other) = delete;
  Flex& operator=(const Flex& other) = delete;
  Flex(Flex&& other) noexcept = default;
  Flex& operator=(Flex&& other) noexcept = default;

  void NextLine();

  void Shrink();

  void Write32(uint32_t w);
  void Write64(uint64_t w);

  void MaskWrite(int count, bool white);

  void MaskFinalise();

  void IncrementLast32(int inc) const;

  uint8_t ReadBackwards8();
  uint16_t ReadBackwards16();
  uint32_t ReadBackwards32();
  uint64_t ReadBackwards64();

  uint32_t ReadForwards32();

  void Copy(uint8_t* src, int len);

  void Start();

  void End();

  std::unique_ptr<uint8_t, FreeDeleter> data_ = nullptr;
  int width_;
  int height_;
  std::vector<int> rows_;
  int y_ = -1;

 private:
  int size_ = 0;
  int p_ = 0;
  int end_p_ = 0;
  int mask_count_ = 0;
  bool mask_white_;
  bool first_;
};

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
