#pragma once

#include <immintrin.h>

namespace multiblend::memory {

// Would use unique_ptr, but then gcc gives warnings:
// "warning: ignoring attributes on template argument"

class AlignedM128Ptr {
  // RAII class for an array of __m128
 public:
  AlignedM128Ptr() : ptr_(nullptr) {}
  AlignedM128Ptr(__m128* ptr) : ptr_(ptr) {}
  ~AlignedM128Ptr();

  AlignedM128Ptr(const AlignedM128Ptr& other) = delete;
  AlignedM128Ptr& operator=(const AlignedM128Ptr& other) = delete;
  AlignedM128Ptr(AlignedM128Ptr&& other) noexcept;
  AlignedM128Ptr& operator=(AlignedM128Ptr&& other) noexcept;

  __m128* get() const { return ptr_; }

  __m128& operator[](int i) { return ptr_[i]; }
  const __m128& operator[](int i) const { return ptr_[i]; }

 private:
  __m128* ptr_ = nullptr;
};
}  // namespace multiblend::memory
