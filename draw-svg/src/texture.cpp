#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Advanced Task
  // Implement mipmap for trilinear filtering

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    /* Advanced Task: Implement MipMap
    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
    */

    MipLevel& mip = tex.mipmap[i];
    MipLevel& mipPrev = tex.mipmap[i-1];
    for(size_t i = 0; i < mip.width; i += 1) {
      for(size_t j=0; j< mip.height; j++){
        Color c;
        Color c1, c2, c3, c4;
        c1.r = mipPrev.texels[4*(2*i + 2*j*mipPrev.width)]/255.0;
        c1.g = mipPrev.texels[4*(2*i + 2*j*mipPrev.width)+1]/255.0;
        c1.b = mipPrev.texels[4*(2*i + 2*j*mipPrev.width)+2]/255.0;
        c1.a = mipPrev.texels[4*(2*i + 2*j*mipPrev.width)+3]/255.0;

        c2.r = mipPrev.texels[4*(2*i + 2*j*mipPrev.width +1)]/255.0;
        c2.g = mipPrev.texels[4*(2*i + 2*j*mipPrev.width +1)+1]/255.0;
        c2.b = mipPrev.texels[4*(2*i + 2*j*mipPrev.width +1)+2]/255.0;
        c2.a = mipPrev.texels[4*(2*i + 2*j*mipPrev.width +1)+3]/255.0;

        c3.r = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width)]/255.0;
        c3.g = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width)+1]/255.0;
        c3.b = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width)+2]/255.0;
        c3.a = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width)+3]/255.0;

        c4.r = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width +1)]/255.0;
        c4.g = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width +1)+1]/255.0;
        c4.b = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width +1)+2]/255.0;
        c4.a = mipPrev.texels[4*(2*i + (2*j+1)*mipPrev.width +1)+3]/255.0;

        c.r = 0.25*c1.r + 0.25*c2.r + 0.25*c3.r + 0.25*c4.r;
        c.g = 0.25*c1.g + 0.25*c2.g + 0.25*c3.g + 0.25*c4.g;
        c.b = 0.25*c1.b + 0.25*c2.b + 0.25*c3.b + 0.25*c4.b;
        c.a = 0.25*c1.a + 0.25*c2.a + 0.25*c3.a + 0.25*c4.a;

        float_to_uint8( &mip.texels[4*(i + j*mip.width)], &c.r);
      }
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {
  // return magenta for invalid level
  if(level<0 || level>=tex.mipmap.size()){
    return Color(1,0,1,1);   
  }
                                  
  // Task 4: Implement nearest neighbour interpolation
  int u_coord = int(u*tex.width), v_coord = int(v*tex.height);
  MipLevel ml = tex.mipmap[level];
  assert(ml.height == tex.height && ml.width == tex.width);
  Color col;
  col.r = ml.texels[4*(u_coord+v_coord*tex.width)]/255.0;
  col.g = ml.texels[4*(u_coord+v_coord*tex.width)+1]/255.0;
  col.b = ml.texels[4*(u_coord+v_coord*tex.width)+2]/255.0;
  col.a = ml.texels[4*(u_coord+v_coord*tex.width)+3]/255.0;
  return col;
}

Color lerp(float x, Color col1, Color col2){
  Color col = (1-x)*col1 + x*col2;
  // col.r = (1-x)*col1.r + x*col2.r;
  // col.g = (1-x)*col1.g + x*col2.g;
  // col.b = (1-x)*col1.b + x*col2.b;
  // col.a = (1-x)*col1.a + x*col2.a; 
  return col;
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
                                      
  
  // Task 4: Implement bilinear filtering
  if(level<0 || level>=tex.mipmap.size()){
    return Color(1,0,1,1);   
  }
  
  MipLevel ml = tex.mipmap[level];

  if (u < 0.0) {
    u = 0.0;
  }
  if (v < 0.0) {
    v = 0.0;
  }
  if (u > 1.0) {
    u = 1.0;
  }
  if (v > 1.0) {
    v = 1.0;
  }
  u*=tex.width;
  v*=tex.height;

  int ul = (u - int(u) >= 0.5)? int(int(u)+0.5) : int(int(u)-0.5), ur = min(ul+1, int(tex.width-1));
  int vt = (v - int(v) >= 0.5)? int(int(v)+0.5) : int(int(v)-0.5), vb = min(vt+1, int(tex.height-1));

  
  // assert(ml.height == tex.height && ml.width == tex.width);

  Color col_lt, col_lb, col_rt, col_rb;
  int idx = 4*(ul+vt*tex.width);
  col_lt.r = ml.texels[idx]/255.0;
  col_lt.g = ml.texels[idx+1]/255.0;
  col_lt.b = ml.texels[idx+2]/255.0;
  col_lt.a = ml.texels[idx+3]/255.0;

  idx = 4*(ul+vb*tex.width);
  col_lb.r = ml.texels[idx]/255.0;
  col_lb.g = ml.texels[idx+1]/255.0;
  col_lb.b = ml.texels[idx+2]/255.0;
  col_lb.a = ml.texels[idx+3]/255.0;

  idx = 4*(ur+vt*tex.width);
  col_rt.r = ml.texels[idx]/255.0;
  col_rt.g = ml.texels[idx+1]/255.0;
  col_rt.b = ml.texels[idx+2]/255.0;
  col_rt.a = ml.texels[idx+3]/255.0;

  idx = 4*(ur+vb*tex.width);
  col_rb.r = ml.texels[idx]/255.0;
  col_rb.g = ml.texels[idx+1]/255.0;
  col_rb.b = ml.texels[idx+2]/255.0;
  col_rb.a = ml.texels[idx+3]/255.0;

  // float s = (u - (ul+0.5))/tex.width, t = (v - (vt+0.5))/tex.height;
  float s = (u - (ul+0.5)), t = (v - (vt+0.5));
  return lerp(t, lerp(s, col_lt, col_rt), lerp(s, col_lb, col_rb));

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Advanced Task
  // Implement trilinear filtering

  float L = max(u_scale, v_scale);
  float d = log2f(L);
  d = max(0.f,d);
  L = max(1.f, L);
  Color mip_low = sample_bilinear(tex,u,v,int(d));
  Color mip_high = sample_bilinear(tex,u,v,int(d)+1);
  
  /* Implementation details: Interpolation of MIP levels can be done 
  either in exponential scale or linear scale */
  float s = (d - int(d)); // Linear scale interpolation
  // float s = (L / (1 << int(d))) -1.0; // Exponential Scale interpolation
  return lerp(s, mip_low, mip_high);
}

} // namespace CS248
