/*
    Bayer pattern debayering using super-pixel method
    Based on QuickFits implementation by Siyu Zhang
*/

#include "debayer.h"

// RGGB pattern:
// R  G
// G  B
void super_pixel_rggb_u16(const uint16_t* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            // R from top-left, G averaged from top-right and bottom-left, B from bottom-right
            newbuf[idx] = (float)buf[cur];
            newbuf[idx + out_plane_size] = ((float)buf[right] + (float)buf[down]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = (float)buf[down_right];
        }
    }
}

void super_pixel_rggb_f32(const float* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            newbuf[idx] = buf[cur];
            newbuf[idx + out_plane_size] = (buf[right] + buf[down]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = buf[down_right];
        }
    }
}

// BGGR pattern:
// B  G
// G  R
void super_pixel_bggr_u16(const uint16_t* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            // R from bottom-right, G averaged, B from top-left
            newbuf[idx] = (float)buf[down_right];
            newbuf[idx + out_plane_size] = ((float)buf[right] + (float)buf[down]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = (float)buf[cur];
        }
    }
}

void super_pixel_bggr_f32(const float* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            newbuf[idx] = buf[down_right];
            newbuf[idx + out_plane_size] = (buf[right] + buf[down]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = buf[cur];
        }
    }
}

// GBRG pattern:
// G  B
// R  G
void super_pixel_gbrg_u16(const uint16_t* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            // R from bottom-left, G averaged, B from top-right
            newbuf[idx] = (float)buf[down];
            newbuf[idx + out_plane_size] = ((float)buf[cur] + (float)buf[down_right]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = (float)buf[right];
        }
    }
}

void super_pixel_gbrg_f32(const float* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            newbuf[idx] = buf[down];
            newbuf[idx + out_plane_size] = (buf[cur] + buf[down_right]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = buf[right];
        }
    }
}

// GRBG pattern:
// G  R
// B  G
void super_pixel_grbg_u16(const uint16_t* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            // R from top-right, G averaged, B from bottom-left
            newbuf[idx] = (float)buf[right];
            newbuf[idx + out_plane_size] = ((float)buf[cur] + (float)buf[down_right]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = (float)buf[down];
        }
    }
}

void super_pixel_grbg_f32(const float* buf, float* newbuf, int width, int height) {
    int out_width = width / 2;
    int out_height = height / 2;
    int out_plane_size = out_width * out_height;

    for (int row = 0, iout = 0; row + 1 < height; row += 2, iout++) {
        for (int col = 0, jout = 0; col + 1 < width; col += 2, jout++) {
            int idx = iout * out_width + jout;
            int cur = row * width + col;
            int right = cur + 1;
            int down = cur + width;
            int down_right = down + 1;

            newbuf[idx] = buf[right];
            newbuf[idx + out_plane_size] = (buf[cur] + buf[down_right]) / 2.0f;
            newbuf[idx + out_plane_size * 2] = buf[down];
        }
    }
}
