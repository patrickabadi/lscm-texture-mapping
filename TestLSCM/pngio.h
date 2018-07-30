#pragma once

#include <stdlib.h>
#include <png.h>

namespace common
{

  typedef enum png_color_type : png_byte {
    GRAY = PNG_COLOR_TYPE_GRAY,
    GRAY_A = PNG_COLOR_TYPE_GRAY_ALPHA,
    PALETTE = PNG_COLOR_TYPE_PALETTE,
    RGB = PNG_COLOR_TYPE_RGB,
    RGB_A = PNG_COLOR_TYPE_RGB_ALPHA
  } png_color_type;

  class pngio
  {
  public:
    pngio (const png_uint_16 width, const png_uint_16 height, png_color_type color_type);
    ~pngio ();

    bool Save (const char const* filename);
    void WriteAt (const png_uint_16 x, const png_uint_16 y, const unsigned char r, const unsigned char g, const unsigned char b);
    void WriteBlockAt(const png_uint_16 x, const png_uint_16 y, int width, int height, unsigned char* data);

  private:
    png_structp png_;
    png_infop png_info_;

    int width_;
    int height_;

    png_byte color_type_;
    png_byte bit_depth_;
    png_bytep *row_pointers_;

    void allocate_memory ();
    void release_memory ();
  };
}