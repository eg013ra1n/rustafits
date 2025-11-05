/*
    Bayer pattern debayering using super-pixel method
    Based on QuickFits implementation by Siyu Zhang
    Adapted for FITS processing
*/

#ifndef DEBAYER_H
#define DEBAYER_H

#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

// Super-pixel debayer for RGGB pattern
// Input: buf (width x height), Output: newbuf (3 planes of width/2 x height/2)
void super_pixel_rggb_u16(const uint16_t* buf, float* newbuf, int width, int height);
void super_pixel_rggb_f32(const float* buf, float* newbuf, int width, int height);

// Super-pixel debayer for BGGR pattern
void super_pixel_bggr_u16(const uint16_t* buf, float* newbuf, int width, int height);
void super_pixel_bggr_f32(const float* buf, float* newbuf, int width, int height);

// Super-pixel debayer for GBRG pattern
void super_pixel_gbrg_u16(const uint16_t* buf, float* newbuf, int width, int height);
void super_pixel_gbrg_f32(const float* buf, float* newbuf, int width, int height);

// Super-pixel debayer for GRBG pattern
void super_pixel_grbg_u16(const uint16_t* buf, float* newbuf, int width, int height);
void super_pixel_grbg_f32(const float* buf, float* newbuf, int width, int height);

#ifdef __cplusplus
}
#endif

#endif /* DEBAYER_H */
