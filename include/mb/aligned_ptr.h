// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cmath>
#include <new>

#include <simde/x86/sse4.1.h>

namespace multiblend::memory {

// Would use unique_ptr, but then gcc gives warnings:
// "warning: ignoring attributes on template argument"

class AlignedM128Ptr {
  // RAII class for an array of simde__m128
 public:
  static constexpr std::align_val_t kAlignment = std::align_val_t{16};

  AlignedM128Ptr() = default;
  ~AlignedM128Ptr();

  AlignedM128Ptr(const AlignedM128Ptr& other) = delete;
  AlignedM128Ptr& operator=(const AlignedM128Ptr& other) = delete;
  AlignedM128Ptr(AlignedM128Ptr&& other) noexcept;
  AlignedM128Ptr& operator=(AlignedM128Ptr&& other) noexcept;

  [[nodiscard]] simde__m128* get() const { return ptr_; }

  simde__m128& operator[](int i) { return ptr_[i]; }
  const simde__m128& operator[](int i) const { return ptr_[i]; }

 private:
  explicit AlignedM128Ptr(simde__m128* ptr) : ptr_(ptr) {}
  simde__m128* ptr_ = nullptr;

  friend AlignedM128Ptr AllocAlignedM128(std::size_t size_bytes);
};

AlignedM128Ptr AllocAlignedM128(std::size_t size_bytes);

}  // namespace multiblend::memory
