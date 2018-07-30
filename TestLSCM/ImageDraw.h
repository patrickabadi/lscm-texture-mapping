#pragma once

#include "pngio.h"
#include "LSCM.h"

using namespace TEX_MAPPER;

class ImageDraw
{
public:
  ImageDraw (common::pngio& png);
  ~ImageDraw ();

  void SetPixel ( const vec2& coord, const color& clr );
  void SetPixel ( png_uint_16 x, png_uint_16 y, const color& clr );
  void SetPixel ( png_uint_16 x, png_uint_16 y, png_byte r, png_byte g, png_byte b );

  void InterpolatedTriangle ( 
    const vec2& v0, const color& c0,
    const vec2& v1, const color& c1,
    const vec2& v2, const color& c2 );
private:
  common::pngio* _png;
};

