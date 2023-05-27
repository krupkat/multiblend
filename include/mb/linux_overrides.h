#pragma once

#ifdef __APPLE__
#define memalign(a, b) malloc((b))
#else
#include <malloc.h>  // memalign
#endif

#ifndef _WIN32

#include <cstring>    // memset
#include <strings.h>  // strcasecmp

inline int _stricmp(const char* a, const char* b) { return strcasecmp(a, b); }

#define ZeroMemory(a, b) memset(a, 0, b)

#define sprintf_s sprintf

#define sscanf_s sscanf

#define strcpy_s(a, b) strcpy(a, b)

inline void* _aligned_malloc(std::size_t size, int boundary) {
  return memalign(boundary, size);
}

inline void _aligned_free(void* a) { free(a); }

inline void fopen_s(FILE** f, const char* filename, const char* mode) {
  *f = fopen(filename, mode);
}

#endif
