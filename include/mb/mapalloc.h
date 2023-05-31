#pragma once

#include <cstddef>
#include <memory>
#include <new>

namespace multiblend::memory {

template <typename TType>
class AlignedAllocDeleter {
 public:
  AlignedAllocDeleter() = default;
  explicit AlignedAllocDeleter(std::align_val_t alignment)
      : alignment_{alignment} {}
  void operator()(TType* ptr) const { ::operator delete(ptr, alignment_); }

 private:
  std::align_val_t alignment_{16};
};

template <typename TType>
using AlignedAllocPtr = std::unique_ptr<TType, AlignedAllocDeleter<TType>>;

template <typename TType>
AlignedAllocPtr<TType> AllocAligned(std::size_t size_bytes,
                                    std::size_t alignment = 16) {
  auto alignment_val = std::align_val_t{alignment};
  return {static_cast<TType*>(::operator new(size_bytes, alignment_val)),
          AlignedAllocDeleter<TType>{alignment_val}};
}

void CacheThreshold(std::size_t limit);
void SetTmpdir(const char* _tmpdir);

// TODO(krupkat): cleanup the original MapAlloc class

/*
class MapAlloc {
 private:
  class MapAllocObject {
   public:
    MapAllocObject(std::size_t size, int alignment);
    ~MapAllocObject();
    void* GetPointer();
    [[nodiscard]] std::size_t GetSize() const { return size_; }
    [[nodiscard]] bool IsFile() const;

   private:
#ifdef _WIN32
    HANDLE file_ = nullptr;
    HANDLE map_ = nullptr;
#else
    int file_ = 0;
#endif
    void* pointer_ = nullptr;
    std::size_t size_;
  };

  static std::vector<MapAllocObject*> objects_;
  static char tmpdir_[256];
  static char filename_[512];
  static int suffix_;
  static std::size_t cache_threshold_;
  static std::size_t total_allocated_;

 public:
  static void* Alloc(std::size_t size, int alignment = 16);
  static void Free(void* p);
  static std::size_t GetSize(void* p);
  static void CacheThreshold(std::size_t limit);
  static void SetTmpdir(const char* _tmpdir);
  static bool LastFile() { return objects_.back()->IsFile(); }
  // static bool last_mapped;
};
*/

}  // namespace multiblend::memory
