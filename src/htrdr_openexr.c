/* Copyright (C) 2018 Université Paul Sabatier, |Meso|Star>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>. */

#include "htrdr_openexr.h"

#include <OpenEXR/Iex.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfIO.h>
#include <OpenEXR/ImfOutputFile.h>

/*******************************************************************************
 * Helper class
 ******************************************************************************/
class OFStream : public Imf::OStream {
public:
  OFStream(FILE* stream, const char* filename)
    : OStream(filename), _stream(stream) {}

  void write(const char c[], int n)
  {
    if(n < 0 || fwrite(c, 1, n, _stream) != (size_t)n) {
      std::stringstream s;
      s << "Error writing the OpenEXR image into`" << fileName() << "'.";

      throw Iex::IoExc(s);
    }
  }

  Imf::Int64 tellp()
  {
    const long l = ftell(_stream);
    if(l == -1) Iex::throwErrnoExc();
    return Imf::Int64(l);
  }

  void seekp(Imf::Int64 pos)
  {
    const int err = fseek(_stream, long(pos), SEEK_SET);
    if(err == -1) Iex::throwErrnoExc();
  }

private:
  FILE* _stream;
};

/*******************************************************************************
 * Local functions
 ******************************************************************************/
res_T
htrdr_openexr_write
  (const float* pixels, 
   const size_t width,
   const size_t height,
   const size_t pitch,
   const char* output_name,
   FILE* output)
{
  res_T res = RES_OK;

  /* Init the OpenEXR slice */
  Imf::Slice slice
    (Imf::FLOAT, (char*)pixels, sizeof(*pixels), pitch);

  /* Setup the OpenEXR framebuffer */
  Imf::FrameBuffer framebuffer;
  framebuffer.insert("L", slice);

  /* Define the OpenEXR header */
  Imf::Header header((int)width, (int)height);
  header.channels().insert("L", Imf::FLOAT);

  /* Setup the OpenEXR output file from the submitted C stream */
  OFStream out_stream(output, output_name);
  Imf::OutputFile out_file(out_stream, header);

  /* Write data */
  out_file.setFrameBuffer(framebuffer);
  try {
    out_file.writePixels((int)height);
  } catch(...) {
    res = RES_IO_ERR;
  }
  return res;
}

