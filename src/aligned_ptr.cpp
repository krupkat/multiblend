#include "mb/aligned_ptr.h"

#include <utility>

#include "mb/linux_overrides.h"

namespace multiblend::memory {

AlignedM128Ptr AllocAlignedM128(std::size_t size_bytes) {
  return AlignedM128Ptr(static_cast<__m128*>(
      ::operator new(size_bytes, AlignedM128Ptr::kAlignment)));
}

AlignedM128Ptr::~AlignedM128Ptr() {
  ::operator delete(ptr_, AlignedM128Ptr::kAlignment);
}

AlignedM128Ptr::AlignedM128Ptr(AlignedM128Ptr&& other) noexcept {
  *this = std::move(other);
}

AlignedM128Ptr& AlignedM128Ptr::operator=(AlignedM128Ptr&& other) noexcept {
  if (this != &other) {
    ::operator delete(ptr_, AlignedM128Ptr::kAlignment);

    ptr_ = other.ptr_;
    other.ptr_ = nullptr;
  }
  return *this;
}

}  // namespace multiblend::memory
