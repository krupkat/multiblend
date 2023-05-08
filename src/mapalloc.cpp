#include "src/mapalloc.h"

#include <cstring>

#ifndef _WIN32
#include <cstdlib>
#include <errno.h>
#include <unistd.h>

#include <sys/mman.h>
#endif

#include "src/linux_overrides.h"

namespace multiblend::memory {

std::vector<MapAlloc::MapAllocObject*> MapAlloc::objects_;
char MapAlloc::tmpdir_[256] = "";
char MapAlloc::filename_[512];
int MapAlloc::suffix_ = 0;
size_t MapAlloc::cache_threshold_ = ~(size_t)0;
size_t MapAlloc::total_allocated_ = 0;

/***********************************************************************
 * MapAlloc
 ***********************************************************************/
void MapAlloc::CacheThreshold(size_t limit) { cache_threshold_ = limit; }

void* MapAlloc::Alloc(size_t size, int alignment) {
  MapAllocObject* m = new MapAllocObject(size, alignment);
  objects_.push_back(m);
  return m->GetPointer();
}

void MapAlloc::Free(void* p) {
  for (auto it = objects_.begin(); it < objects_.end(); ++it) {
    if ((*it)->GetPointer() == p) {
      delete (*it);
      objects_.erase(it);
      break;
    }
  }
}

size_t MapAlloc::GetSize(void* p) {
  for (auto it = objects_.begin(); it < objects_.end(); ++it) {
    if ((*it)->GetPointer() == p) {
      return (*it)->GetSize();
    }
  }

  return 0;
}

void MapAlloc::SetTmpdir(const char* _tmpdir) {
  strcpy_s(tmpdir_, _tmpdir);
  size_t l = strlen(tmpdir_);
  while (tmpdir_[l - 1] == '\\' || tmpdir_[l - 1] == '/' && l > 0)
    tmpdir_[--l] = 0;
}

/***********************************************************************
 * MapAllocObject
 ***********************************************************************/
MapAlloc::MapAllocObject::MapAllocObject(size_t size, int alignment)
    : size_(size) {
  if (total_allocated_ + size_ < cache_threshold_) {
    pointer_ = _aligned_malloc(size_, alignment);
  }

  if (pointer_ == nullptr) {
#ifdef _WIN32
    if (!tmpdir_[0]) {
      GetTempPath(256, tmpdir_);
      size_t l = strlen(tmpdir_);
      while (tmpdir_[l - 1] == '\\' || tmpdir_[l - 1] == '/' && l > 0) {
        tmpdir_[--l] = 0;
      }
    }

    while (true) {
      sprintf_s(filename_, "%s\\_mb%05d.tmp", tmpdir_, suffix_++);
      file_ = CreateFile(filename_, GENERIC_ALL, 0, NULL, CREATE_NEW,
                         FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE |
                             FILE_FLAG_SEQUENTIAL_SCAN,
                         NULL);
      if (file_ != INVALID_HANDLE_VALUE) break;
      if (GetLastError() != 80) {
        char buf[256];
        FormatMessage(
            FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL,
            GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), buf,
            sizeof(buf), NULL);
        sprintf_s(filename_, "Could not create temp file in %s\\: %s", tmpdir_,
                  buf);
        throw(filename_);
      }
      if (suffix_ == 65536) {
        sprintf_s(filename_,
                  "Could not create temp file in %s\\: suffixes exhausted",
                  tmpdir_);
        throw(filename_);
      }
    }

    map_ = CreateFileMapping(file_, NULL, PAGE_READWRITE, size >> 32,
                             size & 0xffffffff, NULL);
    if (!map_) {
      sprintf_s(filename_, "Could not allocate %zu temporary bytes in %s", size,
                tmpdir_);
      throw(filename_);
    }

    pointer_ = MapViewOfFile(map_, FILE_MAP_ALL_ACCESS, 0, 0, 0);
    if (!pointer_) {
      sprintf_s(filename_, "Could not map view of temporary file");
      throw(filename_);
    }
#else
    if (tmpdir_[0] == 0) {
      char* td = getenv("TMPDIR");
      if (td != nullptr) {
        strcpy(tmpdir_, td);
      } else {
        strcpy(tmpdir_, "/tmp");
      }
    }

    sprintf(filename_, "%s/.mbXXXXXX", tmpdir_);
    file_ = mkstemp(filename_);

    if (file_ <= 0) {
      sprintf(filename_, "Could not create temp file in %s/: %s", tmpdir_,
              strerror(errno));
      throw(filename_);
    }

    if (ftruncate(file_, size_) != 0) {
      unlink(filename_);
      sprintf(filename_, "Could not allocate %zu temporary bytes in %s: %s",
              size_, tmpdir_, strerror(errno));
      throw(filename_);
    }

    pointer_ = mmap(NULL, size_, PROT_READ | PROT_WRITE, MAP_PRIVATE, file_, 0);
    if (pointer_ == MAP_FAILED) {
      unlink(filename_);
      pointer_ = NULL;
      sprintf(filename_, "Could not mmap temporary file");
      throw(filename_);
    }

    unlink(filename_);
#endif
  } else {
    total_allocated_ += size_;
  }
}

MapAlloc::MapAllocObject::~MapAllocObject() {
#ifdef _WIN32
  if (!file_) {
    _aligned_free(pointer_);
    total_allocated_ -= size_;
  } else {
    UnmapViewOfFile(pointer_);
    CloseHandle(map_);
    CloseHandle(file_);
  }
#else
  if (file_ == 0) {
    free(pointer_);
  } else {
    munmap(pointer_, size_);
    close(file_);
  }
#endif
}

void* MapAlloc::MapAllocObject::GetPointer() { return pointer_; }

bool MapAlloc::MapAllocObject::IsFile() { return !(file_ == 0); }

}  // namespace multiblend::memory
