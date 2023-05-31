#pragma once

#ifndef _WIN32

#include <cstdio>     // fopen
#include <strings.h>  // strcasecmp

inline int _stricmp(const char* a, const char* b) { return strcasecmp(a, b); }

#define sscanf_s sscanf

inline void fopen_s(FILE** f, const char* filename, const char* mode) {
  *f = fopen(filename, mode);
}

#endif
