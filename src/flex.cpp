// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#include "mb/flex.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>

namespace multiblend::utils {

namespace {
std::unique_ptr<uint8_t, FreeDeleter> SafeMalloc(size_t size) {
  if (auto tmp = std::unique_ptr<uint8_t, FreeDeleter>{(uint8_t*)malloc(size),
                                                       FreeDeleter{}};
      tmp) {
    return tmp;
  }
  throw std::runtime_error("out of memory");
}

void SafeRealloc(std::unique_ptr<uint8_t, FreeDeleter>& data, size_t size) {
  if (auto tmp =
          std::unique_ptr<uint8_t, FreeDeleter>{
              (uint8_t*)realloc(data.get(), size), FreeDeleter{}};
      tmp) {
    data.release();  // NOLINT(bugprone-unused-return-value): used in realloc
    data = std::move(tmp);
    return;
  }
  throw std::runtime_error("out of memory");
}
}  // namespace

/***********************************************************************
 * Flexible data class
 ***********************************************************************/
void Flex::NextLine() {
  end_p_ = p_;
  if (y_ < height_ - 1) {
    rows_[++y_] = p_;
  }

  if (p_ + (width_ << 2) > size_) {
    if (y_ == 0) {
      size_ = (std::max)(height_, 16) << 4;  // was << 2
      data_ = SafeMalloc(size_);
    } else if (y_ < height_) {
      int prev_size = size_;
      int new_size1 = (p_ / y_) * height_ + (width_ << 4);
      int new_size2 = size_ << 1;
      size_ = (std::max)(new_size1, new_size2);
      SafeRealloc(data_, size_);
    }
  }

  MaskFinalise();
  first_ = true;
}

void Flex::Shrink() {
  SafeRealloc(data_, p_);
  end_p_ = p_;
}

void Flex::Write32(uint32_t w) {
  *((uint32_t*)&data_.get()[p_]) = w;
  p_ += 4;
}
void Flex::Write64(uint64_t w) {
  *((uint64_t*)&data_.get()[p_]) = w;
  p_ += 8;
}

void Flex::MaskWrite(int count, bool white) {
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

void Flex::MaskFinalise() {
  if (mask_count_ != 0) {
    Write32((static_cast<int>(mask_white_) << 31) | mask_count_);
  }
}

void Flex::IncrementLast32(int inc) const {
  *((uint32_t*)&data_.get()[p_ - 4]) += inc;
}

uint8_t Flex::ReadBackwards8() { return data_.get()[--p_]; }
uint16_t Flex::ReadBackwards16() { return *((uint16_t*)&data_.get()[p_ -= 2]); }
uint32_t Flex::ReadBackwards32() { return *((uint32_t*)&data_.get()[p_ -= 4]); }
uint64_t Flex::ReadBackwards64() { return *((uint64_t*)&data_.get()[p_ -= 8]); }

uint32_t Flex::ReadForwards32() {
  uint32_t out = *((uint32_t*)&data_.get()[p_]);
  p_ += 4;
  return out;
}

void Flex::Copy(uint8_t* src, int len) {
  memcpy(&data_.get()[p_], src, len);
  p_ += len;
}

void Flex::Start() { p_ = 0; }

void Flex::End() { p_ = end_p_; }

}  // namespace multiblend::utils
