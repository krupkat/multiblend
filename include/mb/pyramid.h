// SPDX-FileCopyrightText: 2021 David Horman
// SPDX-FileCopyrightText: 2023 Tomas Krupka
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstring>
#include <memory>
#include <vector>

#include <simde/x86/sse4.1.h>

#include "mb/aligned_ptr.h"
#include "mb/mapalloc.h"
#include "mb/threadpool.h"

namespace multiblend {

class Pyramid {
 public:
  struct Level {
    int id;
    int width;
    int height;
    int pitch;
    int m128_pitch;
    std::size_t bytes;
    std::shared_ptr<float> data;
    int x;
    int y;
    bool x_shift;
    bool y_shift;
    bool upper_x_shift;
    int upper_m128_pitch;
    std::vector<int> bands;
  };

 private:
  std::vector<Level> levels_;
  std::vector<memory::AlignedM128Ptr> lines_;
  std::vector<float> lut_;
  int lut_bits_ = 0;
  bool lut_gamma_ = false;
  int out_max_;
  mt::ThreadpoolPtr threadpool_;

  void set_lut(int bits, bool gamma);
  void CopyInterleavedThread_8bit(uint8_t* src_p, int step, int pitch, int sy,
                                  int ey);
  void CopyInterleavedThread_16bit(uint16_t* src_p, int step, int pitch, int sy,
                                   int ey);
  void CopyPlanarThread_8bit(uint8_t* src_p, int pitch, bool gamma, int sy,
                             int ey);
  void CopyPlanarThread_16bit(uint16_t* src_p, int pitch, bool gamma, int sy,
                              int ey);
  void CopyPlanarThread_32bit(simde__m128* src_p, int pitch, bool gamma, int sy,
                              int ey);
  static void Subsample_Squeeze(simde__m128* in, simde__m128* Out, int m128_pitch_in,
                                int m128_pitch_out, simde__m128* mul);
  static void ShrinkThread(simde__m128* line, simde__m128* hi, simde__m128* lo,
                           int m128_pitch_hi, int m128_pitch_lo,
                           int first_bad_line, int height_odd, int sy, int ey,
                           bool x_shift, bool y_shift);

  static void Squeeze(simde__m128* line, simde__m128* lo, int m128_pitch_lo,
                      int m128_pitch_hi, simde__m128 final_mul,
                      bool x_shift);  // was __forceinline
  static void LaplaceThreadWrapper(Level* upper_level, Level* lower_level,
                                   int sy, int ey);
  static void LaplaceThread(Level* upper_level, Level* lower_level, int sy,
                            int ey, simde__m128* temp1, simde__m128* temp2,
                            simde__m128* temp3);
  static void FuseThread(simde__m128* a, simde__m128* b, simde__m128* m, int m128_pitch,
                         int sy, int ey, bool pre, int black);
  void LaplaceCollapse(int n_levels, bool Collapse);

  std::size_t total_bytes_ = 0;

 public:
  Pyramid(int width, int height, int _levels, int x, int y,
          mt::ThreadpoolPtr threadpool);

  Pyramid(const Pyramid& other) = delete;
  Pyramid& operator=(const Pyramid& other) = delete;
  Pyramid(Pyramid&& other) = default;
  Pyramid& operator=(Pyramid&& other) = default;

  static int DefaultNumLevels(int width, int height) {
    return 8; /* (int)ceil(log2(max(width, height))); */
  };
  void Copy(uint8_t* src_p, int step, int pitch, bool gamma, int bits);
  void Subsample(int sub_w, int sub_h, Pyramid* source);
  void Shrink();
  void Laplace() { LaplaceCollapse((int)levels_.size(), false); };
  void Collapse() { LaplaceCollapse((int)levels_.size(), true); };
  void Collapse(int n_levels) { LaplaceCollapse(n_levels, true); };
  float Average();
  void Add(float add, int levels = 0x7fff);
  void MultiplyAndAdd(float add, float mul, int levels = 0x7fff);
  void MultiplyAddClamp(float add, float mul, int level = 0x7fff);
  void Multiply(int level, float mul);
  void MultplyByPyramid(Pyramid* b);
  void Fuse(Pyramid* b, Pyramid* mask, bool pre, int black);
  void Fuse(Pyramid* b, float weight);
  void Denoise(int level, float power, bool gamma);
  void Blend(Pyramid* b);
  void BlurX(float radius, Pyramid* transpose);
  void BlurXThread(float radius, Pyramid* transpose, int sy, int ey);
  template <typename T>
  void Out(T dst_p, int pitch, bool gamma, bool dither, bool clamp,
           int level = 0, int step = 0, int offset = 0, bool chroma = false);
  int GetNLevels() { return (int)levels_.size(); };
  int GetPitch(int level = 0) { return levels_[level].pitch; };
  float* GetData(int level = 0) { return levels_[level].data.get(); };
  int GetWidth(int level = 0) { return levels_[level].width; };
  int GetHeight(int level = 0) { return levels_[level].height; };
  int GetX(int level = 0) { return levels_[level].x; };
  int GetY(int level = 0) { return levels_[level].y; };
  static void LaplaceLine2(simde__m128* hi, simde__m128* temp1, simde__m128* temp2,
                           int m128_pitch);
  static void LaplaceLine3(simde__m128* hi, simde__m128* temp1, simde__m128* temp2,
                           simde__m128* temp3, int m128_pitch);
  static void LaplaceExpand(simde__m128* hi, simde__m128* lo, int m128_pitch_hi,
                            int m128_pitch_lo);
  static void LaplaceExpandShifted(simde__m128* hi, simde__m128* lo, int m128_pitch_hi,
                                   int m128_pitch_lo);
  [[nodiscard]] std::size_t GetTotalBytes() const { return total_bytes_; }
  [[nodiscard]] std::vector<Level>& GetLevels() { return levels_; };
  [[nodiscard]] Level& GetLevel(int level) { return levels_[level]; };
  [[nodiscard]] const Level& GetLevel(int level) const {
    return levels_[level];
  };
  void Png(const char* filename);
};

simde__m128* GetLine(const Pyramid::Level& level, int y);

void GetExpandedLine(const Pyramid::Level& level, simde__m128* temp, int y);

}  // namespace multiblend
