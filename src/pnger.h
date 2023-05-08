#pragma once

#include <stdint.h>

#include <png.h>

namespace multiblend::io::png {

class Pnger {
 public:
  Pnger(const char* filename, const char* name, int width, int height, int type,
        int bpp = 8, FILE* file = nullptr, int compression = -1);
  ~Pnger();

  bool Ready() { return !(file_ == nullptr); };
  void WriteRows(uint8_t** rows, int num_rows);
  void Write();
  static void Quick(char* filename, uint8_t* data, int width, int height,
                    int pitch, int type);
  uint8_t* line_;

 private:
  png_structp png_ptr_;
  png_infop info_ptr_;
  static png_color* palette_;
  FILE* file_;
  int y_;
  int height_;
};

}  // namespace multiblend::io::png
