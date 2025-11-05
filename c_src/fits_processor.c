#include "fits_processor.h"
#include "fitsio.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Thread-local error buffer
static __thread char error_buffer[512] = {0};

const char* get_last_error(void) {
    return error_buffer;
}

static void set_error(const char* msg) {
    snprintf(error_buffer, sizeof(error_buffer), "%s", msg);
}

// Forward declarations
static int read_fits_metadata(fitsfile* fptr, FitsMetadata* meta);
static int read_fits_data_u16(fitsfile* fptr, const FitsMetadata* meta, uint16_t** data);
static int read_fits_data_f32(fitsfile* fptr, const FitsMetadata* meta, float** data);
static void downscale_u16(const uint16_t* in, uint16_t* out, size_t w, size_t h, int factor);
static void downscale_f32(const float* in, float* out, size_t w, size_t h, int factor);
static void bin_2x2_float(const float* in, float* out, size_t w, size_t h);
static void super_pixel_debayer_u16(const uint16_t* in, float* out_rgb,
                                    size_t w, size_t h, BayerPattern pattern);
static void super_pixel_debayer_f32(const float* in, float* out_rgb,
                                    size_t w, size_t h, BayerPattern pattern);
static void compute_stretch_params(const float* data, size_t len, float max_input,
                                   float* shadows, float* highlights, float* midtones);
static void apply_stretch(const float* in, uint8_t* out, size_t len,
                         float shadows, float highlights, float midtones, float max_input);

// Comparison function for qsort
static int compare_float(const void* a, const void* b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

int process_fits_file(const char* fits_path, const ProcessConfig* config, ProcessedImage* out_image) {
    fitsfile* fptr = NULL;
    int status = 0;

    memset(out_image, 0, sizeof(ProcessedImage));

    // Open FITS file
    fits_open_file(&fptr, fits_path, READONLY, &status);
    if (status) {
        char fits_err[512];
        fits_get_errstatus(status, fits_err);
        snprintf(error_buffer, sizeof(error_buffer), "Failed to open FITS file: %s (status=%d)", fits_err, status);
        return -1;
    }

    // Read metadata
    FitsMetadata meta;
    if (read_fits_metadata(fptr, &meta) != 0) {
        fits_close_file(fptr, &status);
        return -1;
    }

    // For now, only handle mono 16-bit images (we'll expand this)
    if (meta.channels != 1) {
        set_error("Only mono images supported for now");
        fits_close_file(fptr, &status);
        return -1;
    }

    size_t width = meta.width;
    size_t height = meta.height;

    // Read pixel data
    void* pixel_data = NULL;
    float* float_data = NULL;

    if (meta.dtype == DTYPE_UINT16) {
        uint16_t* u16_data = NULL;
        if (read_fits_data_u16(fptr, &meta, &u16_data) != 0) {
            fits_close_file(fptr, &status);
            return -1;
        }
        pixel_data = u16_data;

        // Apply downscaling if needed
        if (config->downscale_factor > 1) {
            size_t new_w = width / config->downscale_factor;
            size_t new_h = height / config->downscale_factor;
            uint16_t* downscaled = malloc(new_w * new_h * sizeof(uint16_t));
            downscale_u16(u16_data, downscaled, width, height, config->downscale_factor);
            free(u16_data);
            u16_data = downscaled;
            pixel_data = u16_data;
            width = new_w;
            height = new_h;
        }

        // Check for Bayer pattern and debayer
        if (config->apply_debayer && meta.bayer_pattern != BAYER_NONE) {
            // Debayering reduces resolution by 2
            size_t out_w = width / 2;
            size_t out_h = height / 2;
            float_data = malloc(out_w * out_h * 3 * sizeof(float));
            super_pixel_debayer_u16(u16_data, float_data, width, height, meta.bayer_pattern);
            free(u16_data);
            width = out_w;
            height = out_h;
            out_image->is_color = 1;
        } else {
            // Convert to float for stretching
            float_data = malloc(width * height * sizeof(float));
            for (size_t i = 0; i < width * height; i++) {
                float_data[i] = (float)u16_data[i];
            }
            free(u16_data);

            // Apply 2x2 binning in preview mode for mono images
            if (config->preview_mode && meta.bayer_pattern == BAYER_NONE) {
                size_t new_w = width / 2;
                size_t new_h = height / 2;
                float* binned = malloc(new_w * new_h * sizeof(float));
                bin_2x2_float(float_data, binned, width, height);
                free(float_data);
                float_data = binned;
                width = new_w;
                height = new_h;
            }

            out_image->is_color = 0;
        }
    } else {
        // Float data
        float* f32_data = NULL;
        if (read_fits_data_f32(fptr, &meta, &f32_data) != 0) {
            fits_close_file(fptr, &status);
            return -1;
        }

        // Apply downscaling if needed
        if (config->downscale_factor > 1) {
            size_t new_w = width / config->downscale_factor;
            size_t new_h = height / config->downscale_factor;
            float* downscaled = malloc(new_w * new_h * sizeof(float));
            downscale_f32(f32_data, downscaled, width, height, config->downscale_factor);
            free(f32_data);
            f32_data = downscaled;
            width = new_w;
            height = new_h;
        }

        // Check for Bayer and debayer
        if (config->apply_debayer && meta.bayer_pattern != BAYER_NONE) {
            size_t out_w = width / 2;
            size_t out_h = height / 2;
            float_data = malloc(out_w * out_h * 3 * sizeof(float));
            super_pixel_debayer_f32(f32_data, float_data, width, height, meta.bayer_pattern);
            free(f32_data);
            width = out_w;
            height = out_h;
            out_image->is_color = 1;
        } else {
            float_data = f32_data;

            // Apply 2x2 binning in preview mode for mono images
            if (config->preview_mode && meta.bayer_pattern == BAYER_NONE) {
                size_t new_w = width / 2;
                size_t new_h = height / 2;
                float* binned = malloc(new_w * new_h * sizeof(float));
                bin_2x2_float(float_data, binned, width, height);
                free(float_data);
                float_data = binned;
                width = new_w;
                height = new_h;
            }

            out_image->is_color = 0;
        }
    }

    fits_close_file(fptr, &status);

    // Apply stretch and convert to 8-bit
    int num_channels = out_image->is_color ? 3 : 1;
    out_image->width = width;
    out_image->height = height;
    out_image->data = malloc(width * height * 3);  // Always RGB output

    if (config->auto_stretch) {
        for (int c = 0; c < num_channels; c++) {
            float shadows, highlights, midtones;
            size_t channel_size = width * height;
            float* channel_data = float_data + (c * channel_size);

            // Use fixed input range like QuickLook (not actual max)
            float max_input = 65536.0f;  // For 16-bit data

            compute_stretch_params(channel_data, channel_size, max_input,
                                 &shadows, &highlights, &midtones);

            // Apply stretch to this channel - using QuickFits formula exactly
            // Precompute constants
            float hs_range_factor = (highlights == shadows) ? 1.0f : 1.0f / (highlights - shadows);
            float native_shadows = shadows * max_input;
            float native_highlights = highlights * max_input;
            float k1 = (midtones - 1.0f) * hs_range_factor * 255.0f / max_input;
            float k2 = ((2.0f * midtones) - 1.0f) * hs_range_factor / max_input;

            if (out_image->is_color) {
                // Write to separate RGB channels
                for (size_t i = 0; i < channel_size; i++) {
                    float input = channel_data[i];
                    float output;

                    if (input < native_shadows) {
                        output = 0.0f;
                    } else if (input >= native_highlights) {
                        output = 255.0f;
                    } else {
                        float input_floored = input - native_shadows;
                        output = (input_floored * k1) / (input_floored * k2 - midtones);
                    }

                    uint8_t val = (uint8_t)fmaxf(0.0f, fminf(255.0f, output));
                    out_image->data[i * 3 + c] = val;
                }
            } else {
                // Grayscale - write same value to R, G, B
                // Optimized: compute stretch once, then replicate to RGB
                for (size_t i = 0; i < channel_size; i++) {
                    float input = channel_data[i];
                    float output;

                    if (input < native_shadows) {
                        output = 0.0f;
                    } else if (input >= native_highlights) {
                        output = 255.0f;
                    } else {
                        float input_floored = input - native_shadows;
                        output = (input_floored * k1) / (input_floored * k2 - midtones);
                    }

                    uint8_t val = (uint8_t)fmaxf(0.0f, fminf(255.0f, output));

                    // Optimized write: set all 3 bytes at once using memset
                    memset(&out_image->data[i * 3], val, 3);
                }
            }
        }
    }

    free(float_data);

    // Handle vertical flip
    if (meta.flip_vertical) {
        uint8_t* flipped = malloc(width * height * 3);
        for (size_t y = 0; y < height; y++) {
            memcpy(flipped + y * width * 3,
                   out_image->data + (height - 1 - y) * width * 3,
                   width * 3);
        }
        free(out_image->data);
        out_image->data = flipped;
    }

    return 0;
}

static int read_fits_metadata(fitsfile* fptr, FitsMetadata* meta) {
    int status = 0;
    int naxis = 0;
    long naxes[3] = {0};
    int bitpix = 0;

    // Get image dimensions
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 3, naxes, &status);
    fits_get_img_type(fptr, &bitpix, &status);

    if (status) {
        set_error("Failed to read FITS header");
        return -1;
    }

    if (naxis == 2) {
        meta->width = naxes[0];
        meta->height = naxes[1];
        meta->channels = 1;
    } else if (naxis == 3) {
        meta->width = naxes[0];
        meta->height = naxes[1];
        meta->channels = naxes[2];
    } else {
        set_error("Unsupported FITS dimensions");
        return -1;
    }

    // Determine data type
    if (bitpix == USHORT_IMG || bitpix == SHORT_IMG) {
        meta->dtype = DTYPE_UINT16;
    } else if (bitpix == FLOAT_IMG) {
        meta->dtype = DTYPE_FLOAT32;
    } else {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Unsupported BITPIX value: %d", bitpix);
        set_error(err_msg);
        return -1;
    }

    // Check for Bayer pattern
    char bayerpat[64] = {0};
    status = 0;
    fits_read_key(fptr, TSTRING, "BAYERPAT", bayerpat, NULL, &status);
    if (status == 0) {
        if (strcmp(bayerpat, "RGGB") == 0) meta->bayer_pattern = BAYER_RGGB;
        else if (strcmp(bayerpat, "BGGR") == 0) meta->bayer_pattern = BAYER_BGGR;
        else if (strcmp(bayerpat, "GBRG") == 0) meta->bayer_pattern = BAYER_GBRG;
        else if (strcmp(bayerpat, "GRBG") == 0) meta->bayer_pattern = BAYER_GRBG;
        else meta->bayer_pattern = BAYER_NONE;
    } else {
        meta->bayer_pattern = BAYER_NONE;
    }

    // Check row order
    char roworder[64] = {0};
    status = 0;
    fits_read_key(fptr, TSTRING, "ROWORDER", roworder, NULL, &status);
    meta->flip_vertical = (status == 0 && strcmp(roworder, "TOP-DOWN") == 0) ? 1 : 0;

    return 0;
}

static int read_fits_data_u16(fitsfile* fptr, const FitsMetadata* meta, uint16_t** data) {
    int status = 0;
    size_t npixels = meta->width * meta->height * meta->channels;
    int bitpix = 0;

    *data = malloc(npixels * sizeof(uint16_t));
    if (!*data) {
        set_error("Out of memory");
        return -1;
    }

    // Get bitpix to determine signed vs unsigned
    fits_get_img_type(fptr, &bitpix, &status);

    long fpixel[3] = {1, 1, 1};

    // Read image data - use fits_read_img which handles BZERO/BSCALE automatically
    int anynul = 0;
    fits_read_img(fptr, TUSHORT, 1, npixels, NULL, *data, &anynul, &status);

    if (status) {
        char fits_err[512];
        fits_get_errstatus(status, fits_err);
        snprintf(error_buffer, sizeof(error_buffer), "Failed to read FITS pixel data (u16): %s (status=%d, bitpix=%d, npixels=%zu)",
                 fits_err, status, bitpix, npixels);
        free(*data);
        *data = NULL;
        return -1;
    }

    return 0;
}

static int read_fits_data_f32(fitsfile* fptr, const FitsMetadata* meta, float** data) {
    int status = 0;
    size_t npixels = meta->width * meta->height * meta->channels;

    *data = malloc(npixels * sizeof(float));
    if (!*data) {
        set_error("Out of memory");
        return -1;
    }

    long fpixel[3] = {1, 1, 1};
    fits_read_pix(fptr, TFLOAT, fpixel, npixels, NULL, *data, NULL, &status);

    if (status) {
        set_error("Failed to read FITS pixel data");
        free(*data);
        *data = NULL;
        return -1;
    }

    return 0;
}

static void downscale_u16(const uint16_t* in, uint16_t* out, size_t w, size_t h, int factor) {
    size_t new_w = w / factor;
    size_t new_h = h / factor;

    for (size_t y = 0; y < new_h; y++) {
        for (size_t x = 0; x < new_w; x++) {
            out[y * new_w + x] = in[(y * factor) * w + (x * factor)];
        }
    }
}

static void downscale_f32(const float* in, float* out, size_t w, size_t h, int factor) {
    size_t new_w = w / factor;
    size_t new_h = h / factor;

    for (size_t y = 0; y < new_h; y++) {
        for (size_t x = 0; x < new_w; x++) {
            out[y * new_w + x] = in[(y * factor) * w + (x * factor)];
        }
    }
}

// 2x2 binning (averaging) for mono images in preview mode
static void bin_2x2_float(const float* in, float* out, size_t w, size_t h) {
    size_t out_w = w / 2;
    size_t out_h = h / 2;

    for (size_t y = 0; y < out_h; y++) {
        for (size_t x = 0; x < out_w; x++) {
            size_t in_y = y * 2;
            size_t in_x = x * 2;

            // Average 2x2 block
            float p00 = in[in_y * w + in_x];
            float p01 = in[in_y * w + in_x + 1];
            float p10 = in[(in_y + 1) * w + in_x];
            float p11 = in[(in_y + 1) * w + in_x + 1];

            out[y * out_w + x] = (p00 + p01 + p10 + p11) * 0.25f;
        }
    }
}

static void super_pixel_debayer_u16(const uint16_t* in, float* out_rgb,
                                    size_t w, size_t h, BayerPattern pattern) {
    size_t out_w = w / 2;
    size_t out_h = h / 2;

    for (size_t y = 0; y < out_h; y++) {
        for (size_t x = 0; x < out_w; x++) {
            size_t in_y = y * 2;
            size_t in_x = x * 2;

            float p00 = (float)in[in_y * w + in_x];
            float p01 = (float)in[in_y * w + in_x + 1];
            float p10 = (float)in[(in_y + 1) * w + in_x];
            float p11 = (float)in[(in_y + 1) * w + in_x + 1];

            float r, g, b;

            switch (pattern) {
                case BAYER_RGGB:
                    r = p00;
                    g = (p01 + p10) / 2.0f;
                    b = p11;
                    break;
                case BAYER_BGGR:
                    b = p00;
                    g = (p01 + p10) / 2.0f;
                    r = p11;
                    break;
                case BAYER_GBRG:
                    b = p01;
                    g = (p00 + p11) / 2.0f;
                    r = p10;
                    break;
                case BAYER_GRBG:
                    r = p01;
                    g = (p00 + p11) / 2.0f;
                    b = p10;
                    break;
                default:
                    r = g = b = 0.0f;
            }

            // Store in planar format: all R, then all G, then all B
            size_t idx = y * out_w + x;
            size_t plane_size = out_w * out_h;
            out_rgb[idx] = r;                       // R plane
            out_rgb[idx + plane_size] = g;          // G plane
            out_rgb[idx + plane_size * 2] = b;      // B plane
        }
    }
}

static void super_pixel_debayer_f32(const float* in, float* out_rgb,
                                    size_t w, size_t h, BayerPattern pattern) {
    size_t out_w = w / 2;
    size_t out_h = h / 2;

    for (size_t y = 0; y < out_h; y++) {
        for (size_t x = 0; x < out_w; x++) {
            size_t in_y = y * 2;
            size_t in_x = x * 2;

            float p00 = in[in_y * w + in_x];
            float p01 = in[in_y * w + in_x + 1];
            float p10 = in[(in_y + 1) * w + in_x];
            float p11 = in[(in_y + 1) * w + in_x + 1];

            float r, g, b;

            switch (pattern) {
                case BAYER_RGGB:
                    r = p00;
                    g = (p01 + p10) / 2.0f;
                    b = p11;
                    break;
                case BAYER_BGGR:
                    b = p00;
                    g = (p01 + p10) / 2.0f;
                    r = p11;
                    break;
                case BAYER_GBRG:
                    b = p01;
                    g = (p00 + p11) / 2.0f;
                    r = p10;
                    break;
                case BAYER_GRBG:
                    r = p01;
                    g = (p00 + p11) / 2.0f;
                    b = p10;
                    break;
                default:
                    r = g = b = 0.0f;
            }

            // Store in planar format: all R, then all G, then all B
            size_t idx = y * out_w + x;
            size_t plane_size = out_w * out_h;
            out_rgb[idx] = r;                       // R plane
            out_rgb[idx + plane_size] = g;          // G plane
            out_rgb[idx + plane_size * 2] = b;      // B plane
        }
    }
}

static void compute_stretch_params(const float* data, size_t len, float max_input,
                                   float* shadows, float* highlights, float* midtones) {
    // Sample up to 500k pixels
    const size_t MAX_SAMPLES = 500000;
    size_t num_samples = (len <= MAX_SAMPLES) ? len : MAX_SAMPLES;

    float* samples = malloc(num_samples * sizeof(float));
    if (len <= MAX_SAMPLES) {
        memcpy(samples, data, num_samples * sizeof(float));
    } else {
        size_t step = len / MAX_SAMPLES;
        for (size_t i = 0; i < num_samples; i++) {
            samples[i] = data[i * step];
        }
    }

    // Compute median
    qsort(samples, num_samples, sizeof(float), compare_float);
    float median = samples[num_samples / 2];

    // Compute MADN (Median Absolute Deviation Normalized)
    float* deviations = malloc(num_samples * sizeof(float));
    for (size_t i = 0; i < num_samples; i++) {
        deviations[i] = fabsf(samples[i] - median);
    }
    qsort(deviations, num_samples, sizeof(float), compare_float);
    float madn = 1.4826f * deviations[num_samples / 2];

    free(samples);
    free(deviations);

    // Normalize
    float norm_median = median / max_input;
    float norm_madn = madn / max_input;

    // Compute shadows and highlights using QuickFits formula (exactly as in Stretch.h)
    bool upper_half = norm_median > 0.5f;

    // Shadows
    if (upper_half || norm_madn == 0.0f) {
        *shadows = 0.0f;
    } else {
        *shadows = fminf(1.0f, fmaxf(0.0f, norm_median + (-2.8f * norm_madn)));
    }

    // Highlights
    if (!upper_half || norm_madn == 0.0f) {
        *highlights = 1.0f;
    } else {
        *highlights = fminf(1.0f, fmaxf(0.0f, norm_median - (-2.8f * norm_madn)));
    }

    // Compute midtones using QuickFits formula (exactly as in Stretch.h)
    const float B = 0.25f;
    float X, M;

    if (!upper_half) {
        X = norm_median - *shadows;
        M = B;
    } else {
        X = B;
        M = *highlights - norm_median;
    }

    // Calculate midtones
    if (X == 0.0f) {
        *midtones = 0.0f;
    } else if (X == M) {
        *midtones = 0.5f;
    } else if (X == 1.0f) {
        *midtones = 1.0f;
    } else {
        *midtones = ((M - 1.0f) * X) / ((2.0f * M - 1.0f) * X - M);
    }

    // Debug output
    fprintf(stderr, "STRETCH DEBUG: median=%.1f, max_input=%.1f, norm_median=%.6f, norm_madn=%.6f\n",
            median, max_input, norm_median, norm_madn);
    fprintf(stderr, "STRETCH DEBUG: shadows=%.6f, highlights=%.6f, midtones=%.6f\n",
            *shadows, *highlights, *midtones);
}

void free_processed_image(ProcessedImage* image) {
    if (image && image->data) {
        free(image->data);
        image->data = NULL;
    }
}
