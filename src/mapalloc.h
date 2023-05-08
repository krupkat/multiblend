#pragma once

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
    MapAllocObject(size_t size, int alignment);
    ~MapAllocObject();
    void* GetPointer();
    [[nodiscard]] size_t GetSize() const { return size_; }
    [[nodiscard]] bool IsFile() const;

   private:
#ifdef _WIN32
    HANDLE file_ = nullptr;
    HANDLE map_ = nullptr;
#else
    int file_ = 0;
#endif
    void* pointer_ = nullptr;
    size_t size_;
  };

  static std::vector<MapAllocObject*> objects_;
  static char tmpdir_[256];
  static char filename_[512];
  static int suffix_;
  static size_t cache_threshold_;
  static size_t total_allocated_;

 public:
  static void* Alloc(size_t size, int alignment = 16);
  static void Free(void* p);
  static size_t GetSize(void* p);
  static void CacheThreshold(size_t limit);
  static void SetTmpdir(const char* _tmpdir);
  static bool LastFile() { return objects_.back()->IsFile(); }
  // static bool last_mapped;
};

}  // namespace multiblend::memory
