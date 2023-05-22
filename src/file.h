#pragma once

#include <cstdio>
#include <memory>

namespace multiblend::io {
class FileDeleter {
 public:
  void operator()(FILE* file) const { fclose(file); }
};

using FilePtr = std::unique_ptr<FILE, FileDeleter>;

}  // namespace multiblend::io
