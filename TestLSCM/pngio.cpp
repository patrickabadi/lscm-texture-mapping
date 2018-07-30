#include <stdlib.h>
#include <png.h>

#include "Helpers.h"
#include "pngio.h"

common::pngio::pngio (const png_uint_16 width, const png_uint_16 height, png_color_type color_type) :
  width_ (width), height_ (height), color_type_ (color_type), png_ (nullptr), png_info_ (nullptr), row_pointers_ (nullptr)
{
  allocate_memory ();
}

common::pngio::~pngio ()
{
  release_memory ();
}

bool common::pngio::Save (const char const * filename)
{
  if (!png_)
  {
    DebugOut ("Save failed: png_ not initialized.");
    return false;
  }

  if (!png_info_)
  {
    DebugOut ("Save failed: png_info_ not initialized.");
    return false;
  }

  FILE *fp = fopen (filename, "wb");
  if (!fp)
  {
    DebugOut ("Failed to open %s for writing", filename);
    return false;
  }

  if (setjmp (png_jmpbuf (png_)))
  {
    png_destroy_write_struct (&png_, &png_info_);
    fclose (fp);
    return false;
  }

  png_set_check_for_invalid_index (png_, 0);

  png_init_io (png_, fp);

  png_write_info(png_, png_info_);

  png_write_image (png_, row_pointers_);
  png_write_end (png_, NULL);

  fclose (fp);

  return true;
}

void common::pngio::WriteAt (const png_uint_16 x, const png_uint_16 y, const unsigned char R, const unsigned char G, const unsigned char B)
{
  if (y >= height_ || x >= width_)
    return;

  png_bytep row = row_pointers_[y];
  png_bytep px = &(row[x * 4]);

  px[0] = (png_byte)R;
  px[1] = (png_byte)G;
  px[2] = (png_byte)B;
}

void common::pngio::WriteBlockAt (const png_uint_16 x, const png_uint_16 y, int blockWidth, int blockHeight, unsigned char* data)
{
  if (y + blockHeight > height_ || x + blockWidth > width_)
  {
    DebugOut ("Requested block write exceeds dimensions of PNG");
    return;
  }

  int idx;
  for (int i = 0; i < blockHeight; i++) 
  {
    png_bytep row = row_pointers_[y+blockHeight-1-i];
    for (int j=0; j< blockWidth; j++ )
    {
      idx = 4 * (i * blockWidth + j);
      // added i*j for debug black lines
      auto r = data[idx];
      auto g = data[idx + 1];
      auto b = data[idx + 2];
      auto a = data[idx + 3];

      png_bytep px = &(row[(x+j) * 4]);
      px[0] = (png_byte)r;
      px[1] = (png_byte)g;
      px[2] = (png_byte)b;
      px[3] = (png_byte)a;
    }
  }
}

void common::pngio::allocate_memory ()
{
  png_ = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_)
    return;

  png_info_ = png_create_info_struct (png_);
  if (!png_info_)
    return;

  png_set_IHDR (
    png_,
    png_info_,
    width_,
    height_,
    8,
    color_type_,
    PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_DEFAULT,
    PNG_FILTER_TYPE_DEFAULT);

  row_pointers_ = (png_bytep*)malloc (sizeof (png_bytep) * height_);

  auto row_bytes = png_get_rowbytes (png_, png_info_);
  for (int y = 0; y < height_; y++) 
  {
    row_pointers_[y] = (png_byte*)malloc (row_bytes);
  }

}

void common::pngio::release_memory ()
{
  if (!row_pointers_)
    return;

  for (int y = 0; y < height_; y++) {
    if(row_pointers_[y])
      free (row_pointers_[y]);
  }
  free (row_pointers_);
}
