#pragma once

#include <memory>
#include <vector>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

namespace multiblend::memory {

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

class MapAllocDeleter {
 public:
  void operator()(void* p) const { MapAlloc::Free(p); }
};

using MapAllocPtr = std::unique_ptr<void, MapAllocDeleter>;

}  // namespace multiblend::memory
