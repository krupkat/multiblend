#pragma once

#include <chrono>
#include <cstring>
#include <stdint.h>
#include <vector>

#include "src/pyramid.h"

namespace multiblend {

extern int verbosity;

/***********************************************************************
 * Flexible data class
 ***********************************************************************/
class Flex {
 public:
  Flex() : width(0), height(0) {}
  Flex(int _width, int _height) : width(_width), height(_height) {
    rows = new int[height];
    NextLine();
  };

  Flex(const Flex& other) = delete;
  Flex& operator=(const Flex& other) = delete;

  Flex(Flex&& other) { *this = std::move(other); }
  Flex& operator=(Flex&& other) {
    if (this != &other) {
      if (data) {
        free(data);
      }
      if (rows) {
        delete[] rows;
      }

      data = other.data;
      width = other.width;
      height = other.height;
      rows = other.rows;
      y = other.y;

      size = other.size;
      p = other.p;
      end_p = other.end_p;
      mask_count = other.mask_count;
      mask_white = other.mask_white;
      first = other.first;

      other.data = nullptr;
      other.rows = nullptr;
    }
    return *this;
  }

  void NextLine() {
    end_p = p;
    if (y < height - 1) rows[++y] = p;

    if (p + (width << 2) > size) {
      if (y == 0) {
        size = (std::max)(height, 16) << 4;  // was << 2
        data = (uint8_t*)malloc(size);
      } else if (y < height) {
        int prev_size = size;
        int new_size1 = (p / y) * height + (width << 4);
        int new_size2 = size << 1;
        size = (std::max)(new_size1, new_size2);
        data = (uint8_t*)realloc(data, size);
      }
    }

    MaskFinalise();
    first = true;
  }

  void Shrink() {
    data = (uint8_t*)realloc(data, p);
    end_p = p;
  }

  void Write32(uint32_t w) {
    *((uint32_t*)&data[p]) = w;
    p += 4;
  }
  void Write64(uint64_t w) {
    *((uint64_t*)&data[p]) = w;
    p += 8;
  }

  void MaskWrite(int count, bool white) {
    if (first) {
      mask_count = count;
      mask_white = white;
      first = false;
    } else {
      if (white == mask_white) {
        mask_count += count;
      } else {
        Write32((mask_white << 31) | mask_count);
        mask_count = count;
        mask_white = white;
      }
    }
  }

  void MaskFinalise() {
    if (mask_count) Write32((mask_white << 31) | mask_count);
  }

  void IncrementLast32(int inc) { *((uint32_t*)&data[p - 4]) += inc; }

  uint8_t ReadBackwards8() { return data[--p]; }
  uint16_t ReadBackwards16() { return *((uint16_t*)&data[p -= 2]); }
  uint32_t ReadBackwards32() { return *((uint32_t*)&data[p -= 4]); }
  uint64_t ReadBackwards64() { return *((uint64_t*)&data[p -= 8]); }

  uint32_t ReadForwards32() {
    uint32_t out = *((uint32_t*)&data[p]);
    p += 4;
    return out;
  }

  void Copy(uint8_t* src, int len) {
    memcpy(&data[p], src, len);
    p += len;
  }

  void Start() { p = 0; }

  void End() { p = end_p; }

  ~Flex() {
    if (data) {
      free(data);
    }
    if (rows) {
      delete[] rows;
    }
  }

  uint8_t* data = NULL;
  int width;
  int height;
  int* rows = nullptr;
  int y = -1;

 private:
  int size = 0;
  int p = 0;
  int end_p = 0;
  int mask_count = 0;
  bool mask_white;
  bool first;
};

/***********************************************************************
 * Timer
 ***********************************************************************/
class Timer {
 public:
  void Start() { start_time = std::chrono::high_resolution_clock::now(); };

  double Read() {
    std::chrono::duration<double> elapsed =
        std::chrono::high_resolution_clock::now() - start_time;
    return elapsed.count();
  };

  void Report(const char* name) { printf("%s: %.3fs\n", name, Read()); };

 private:
  std::chrono::high_resolution_clock::time_point start_time;
};

void Output(int level, const char* fmt, ...);

void die(const char* error, ...);

void ShrinkMasks(std::vector<Flex*>& masks, int n_levels);

void CompositeLine(float* input_p, float* output_p, int i, int x_offset,
                   int in_level_width, int out_level_width, int out_level_pitch,
                   uint8_t* _mask, size_t mask_p);

/***********************************************************************
 * Seam macro
 ***********************************************************************/
// Record macro
#define RECORD(I, C)                                                \
  if ((I) != current_i) {                                           \
    if (mc > 0) {                                                   \
      if (seam_map) memset(&seam_map->line[x - mc], current_i, mc); \
      for (int i = 0; i < n_images; ++i) {                          \
        if (i == current_i) {                                       \
          images[i]->masks[0]->Write32(0xc0000000 | mc);            \
        } else if (i == prev_i || prev_i == -1) {                   \
          images[i]->masks[0]->Write32(0x80000000 | mc);            \
        } else {                                                    \
          images[i]->masks[0]->IncrementLast32(mc);                 \
        }                                                           \
      }                                                             \
    }                                                               \
    prev_i = current_i;                                             \
    mc = C;                                                         \
    current_i = I;                                                  \
  } else {                                                          \
    mc += C;                                                        \
  }
// end macro

void ReadInpaintDT(Flex* flex, int& current_count, int& current_step,
                   uint32_t& dt_val);

void ReadSeamDT(Flex* flex, int& current_count, int64_t& current_step,
                uint64_t& dt_val);

#define SEAM_DT ReadSeamDT(seam_flex, current_count, current_step, dt_val);
#define INPAINT_DT ReadInpaintDT(dt, current_count, current_step, dt_val);

int CompressDTLine(uint32_t* input, uint8_t* output, int width);

int CompressSeamLine(uint64_t* input, uint8_t* output, int width);

void SwapH(Pyramid* py);

void UnswapH(Pyramid* py);

void SwapV(Pyramid* py);

void UnswapV(Pyramid* py);

}  // namespace multiblend
