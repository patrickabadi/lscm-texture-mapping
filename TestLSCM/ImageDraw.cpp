#include "ImageDraw.h"



ImageDraw::ImageDraw (common::pngio& png)
  : _png (&png)
{
}


ImageDraw::~ImageDraw ()
{
}

void ImageDraw::SetPixel ( const vec2& coord, const color& clr )
{
  SetPixel ( (png_uint_16)coord.x, (png_uint_16)coord.y, clr.R, clr.G, clr.B );
}

void ImageDraw::SetPixel ( png_uint_16 x, png_uint_16 y, const color& clr )
{
  SetPixel ( x, y, clr.R, clr.G, clr.B );
}

void ImageDraw::SetPixel ( png_uint_16 x, png_uint_16 y, png_byte r, png_byte g, png_byte b )
{
  _png->WriteAt ( x, y, r, g, b );
}

double edgeFunction ( const vec2 &a, const vec2 &b, const vec2 &c )
{
  return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}

double minX ( const std::vector<vec2> varray, int index = -1, double minx = std::numeric_limits<double>::max () )
{
  return ++index >= varray.size() ? minx : minX ( varray, index, std::min ( minx, varray[index].x ) );
}

double minY ( const std::vector<vec2> varray, int index = -1, double miny = std::numeric_limits<double>::max () )
{
  return ++index >= varray.size () ? miny : minY ( varray, index, std::min ( miny, varray[index].y ) );
}

double maxX ( const std::vector<vec2> varray, int index = -1, double maxx = std::numeric_limits<double>::min () )
{
  return ++index >= varray.size () ? maxx : maxX ( varray, index, std::max ( maxx, varray[index].x ) );
}

double maxY ( const std::vector<vec2> varray, int index = -1, double maxy = std::numeric_limits<double>::min () )
{
  return ++index >= varray.size () ? maxy : maxY ( varray, index, std::max ( maxy, varray[index].y ) );
}

// https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/rasterization-stage

void ImageDraw::InterpolatedTriangle (
  const vec2& v0, const color& c0,
  const vec2& v1, const color& c1,
  const vec2& v2, const color& c2 )
{
  std::vector<vec2> arr = { v0, v1, v2 };

  auto minx = (uint64_t) minX ( arr );
  auto miny = (uint64_t) minY ( arr );
  auto maxx = (uint64_t) maxX ( arr );
  auto maxy = (uint64_t) maxY ( arr );

  auto area = edgeFunction ( v0, v1, v2 );

  for (uint64_t j = miny; j <= maxy; ++j) 
  {
    for (uint64_t i = minx; i <= maxx; ++i) 
    {
      auto p = vec2( i + 0.5f, j + 0.5f );
      auto w0 = edgeFunction ( v1, v2, p );
      auto w1 = edgeFunction ( v2, v0, p );
      auto w2 = edgeFunction ( v0, v1, p );
      if (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0) 
      {
        w0 /= area;
        w1 /= area;
        w2 /= area;
        auto r = w0 * c0.Rf + w1 * c1.Rf + w2 * c2.Rf;
        auto g = w0 * c0.Gf + w1 * c1.Gf + w2 * c2.Gf;
        auto b = w0 * c0.Bf + w1 * c1.Bf + w2 * c2.Bf;
        SetPixel ( i, j, (unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255) );
      }
    }
  }
}
