#pragma once

#include <chrono>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
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

  void NextLine() {
    end_p_ = p_;
    if (y_ < height_ - 1) {
      rows_[++y_] = p_;
    }

    if (p_ + (width_ << 2) > size_) {
      if (y_ == 0) {
        size_ = (std::max)(height_, 16) << 4;  // was << 2
        data_ = {(uint8_t*)malloc(size_), FreeDeleter{}};
      } else if (y_ < height_) {
        int prev_size = size_;
        int new_size1 = (p_ / y_) * height_ + (width_ << 4);
        int new_size2 = size_ << 1;
        size_ = (std::max)(new_size1, new_size2);
        data_ = {(uint8_t*)realloc(data_.release(), size_), FreeDeleter{}};
      }
    }

    MaskFinalise();
    first_ = true;
  }

  void Shrink() {
    data_ = {(uint8_t*)realloc(data_.release(), p_), FreeDeleter{}};
    end_p_ = p_;
  }

  void Write32(uint32_t w) {
    *((uint32_t*)&data_.get()[p_]) = w;
    p_ += 4;
  }
  void Write64(uint64_t w) {
    *((uint64_t*)&data_.get()[p_]) = w;
    p_ += 8;
  }

  void MaskWrite(int count, bool white) {
    if (first_) {
      mask_count_ = count;
      mask_white_ = white;
      first_ = false;
    } else {
      if (white == mask_white_) {
        mask_count_ += count;
      } else {
        Write32((static_cast<int>(mask_white_) << 31) | mask_count_);
        mask_count_ = count;
        mask_white_ = white;
      }
    }
  }

  void MaskFinalise() {
    if (mask_count_ != 0) {
      Write32((static_cast<int>(mask_white_) << 31) | mask_count_);
    }
  }

  void IncrementLast32(int inc) const {
    *((uint32_t*)&data_.get()[p_ - 4]) += inc;
  }

  uint8_t ReadBackwards8() { return data_.get()[--p_]; }
  uint16_t ReadBackwards16() { return *((uint16_t*)&data_.get()[p_ -= 2]); }
  uint32_t ReadBackwards32() { return *((uint32_t*)&data_.get()[p_ -= 4]); }
  uint64_t ReadBackwards64() { return *((uint64_t*)&data_.get()[p_ -= 8]); }

  uint32_t ReadForwards32() {
    uint32_t out = *((uint32_t*)&data_.get()[p_]);
    p_ += 4;
    return out;
  }

  void Copy(uint8_t* src, int len) {
    memcpy(&data_.get()[p_], src, len);
    p_ += len;
  }

  void Start() { p_ = 0; }

  void End() { p_ = end_p_; }

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
