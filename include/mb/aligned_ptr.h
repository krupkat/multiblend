// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <immintrin.h>
#include <new>

namespace multiblend::memory {

// Would use unique_ptr, but then gcc gives warnings:
// "warning: ignoring attributes on template argument"

class AlignedM128Ptr {
  // RAII class for an array of __m128
 public:
  static constexpr std::align_val_t kAlignment = std::align_val_t{16};

  AlignedM128Ptr() = default;
  ~AlignedM128Ptr();

  AlignedM128Ptr(const AlignedM128Ptr& other) = delete;
  AlignedM128Ptr& operator=(const AlignedM128Ptr& other) = delete;
  AlignedM128Ptr(AlignedM128Ptr&& other) noexcept;
  AlignedM128Ptr& operator=(AlignedM128Ptr&& other) noexcept;

  [[nodiscard]] __m128* get() const { return ptr_; }

  __m128& operator[](int i) { return ptr_[i]; }
  const __m128& operator[](int i) const { return ptr_[i]; }

 private:
  explicit AlignedM128Ptr(__m128* ptr) : ptr_(ptr) {}
  __m128* ptr_ = nullptr;

  friend AlignedM128Ptr AllocAlignedM128(std::size_t size_bytes);
};

AlignedM128Ptr AllocAlignedM128(std::size_t size_bytes);

}  // namespace multiblend::memory
