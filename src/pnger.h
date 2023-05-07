#pragma once

#include <stdint.h>

#include <png.h>

namespace multiblend::io::png {

class Pnger {
 public:
  Pnger(const char* filename, const char* name, int w, int _h, int type,
        int bpp = 8, FILE* _f = NULL, int compression = -1);
  ~Pnger();
  bool Ready() { return !!f; };
  void WriteRows(uint8_t** rows, int num_rows);
  void Write();
  static void Quick(char* filename, uint8_t* data, int width, int height,
                    int pitch, int type);
  uint8_t* line;

 private:
  png_structp png_ptr;
  png_infop info_ptr;
  static png_color* palette;
  FILE* f;
  int y;
  int h;
};

}  // namespace multiblend::io::png
