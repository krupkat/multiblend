#include "src/aligned_ptr.h"

#include <utility>

#include "src/linux_overrides.h"

namespace multiblend::memory {

AlignedM128Ptr::~AlignedM128Ptr() {
  if (ptr_) {
    _aligned_free(ptr_);
  }
}

AlignedM128Ptr::AlignedM128Ptr(AlignedM128Ptr&& other) noexcept {
  *this = std::move(other);
}

AlignedM128Ptr& AlignedM128Ptr::operator=(AlignedM128Ptr&& other) noexcept {
  if (this != &other) {
    if (!ptr_) {
      _aligned_free(ptr_);
    }
    ptr_ = other.ptr_;
    other.ptr_ = nullptr;
  }
  return *this;
}

}  // namespace multiblend::memory
