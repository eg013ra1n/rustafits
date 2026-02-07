#include "fits_processor.h"
#include "xisf_reader.h"
#include "fits_reader.h"
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// SIMD support detection
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #include <emmintrin.h>  // SSE2
    #define HAS_SIMD_SSE2 1
#elif defined(__ARM_NEON) || defined(__aarch64__)
    #include <arm_neon.h>
    #define HAS_SIMD_NEON 1
#endif

// Thread-local error buffer (also used by xisf_reader.c)
THREAD_LOCAL char error_buffer[512] = {0};

const char* get_last_error(void) {
    return error_buffer;
}

static void set_error(const char* msg) {
    snprintf(error_buffer, sizeof(error_buffer), "%s", msg);
}

// Forward declarations
static void downscale_u16(const uint16_t* in, uint16_t* out, size_t w, size_t h, int factor);
static void downscale_f32(const float* in, float* out, size_t w, size_t h, int factor);
static void bin_2x2_float(const float* in, float* out, size_t w, size_t h);
static void super_pixel_debayer_u16(const uint16_t* in, float* out_rgb,
                                    size_t w, size_t h, BayerPattern pattern);
static void super_pixel_debayer_f32(const float* in, float* out_rgb,
                                    size_t w, size_t h, BayerPattern pattern);
static void compute_stretch_params(const float* data, size_t len, float max_input,
                                   float* shadows, float* highlights, float* midtones);

// Quickselect algorithm for finding k-th element (faster than full sort for median)
static float quickselect(float* arr, size_t left, size_t right, size_t k) {
    while (left < right) {
        // Use median-of-three for pivot selection
        size_t mid = left + (right - left) / 2;

        // Sort left, mid, right
        if (arr[mid] < arr[left]) {
            float tmp = arr[left]; arr[left] = arr[mid]; arr[mid] = tmp;
        }
        if (arr[right] < arr[left]) {
            float tmp = arr[left]; arr[left] = arr[right]; arr[right] = tmp;
        }
        if (arr[right] < arr[mid]) {
            float tmp = arr[mid]; arr[mid] = arr[right]; arr[right] = tmp;
        }

        // Use middle element as pivot
        float pivot = arr[mid];

        // Move pivot to end
        arr[mid] = arr[right - 1];
        arr[right - 1] = pivot;

        // Partition
        size_t i = left;
        size_t j = right - 1;

        while (1) {
            while (arr[++i] < pivot);
            while (arr[--j] > pivot);
            if (i >= j) break;
            float tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }

        // Restore pivot
        arr[right - 1] = arr[i];
        arr[i] = pivot;

        // Recurse on appropriate partition
        if (i == k) {
            return arr[k];
        } else if (i > k) {
            right = i - 1;
        } else {
            left = i + 1;
        }
    }
    return arr[k];
}

// Find median using quickselect (O(n) average case vs O(n log n) for qsort)
static float find_median(float* data, size_t len) {
    return quickselect(data, 0, len - 1, len / 2);
}

// SIMD-optimized grayscale to RGB replication
static void replicate_gray_to_rgb(const uint8_t* gray, uint8_t* rgb, size_t count) {
    size_t i = 0;
#ifdef HAS_SIMD_SSE2
    // Process 4 pixels at a time: 4 gray bytes -> 12 RGB bytes
    for (; i + 3 < count; i += 4) {
        uint8_t g0 = gray[i], g1 = gray[i+1], g2 = gray[i+2], g3 = gray[i+3];
        uint8_t* dst = &rgb[i * 3];
        // Write 12 bytes: R0G0B0 R1G1B1 R2G2B2 R3G3B3
        dst[0] = g0; dst[1] = g0; dst[2] = g0;
        dst[3] = g1; dst[4] = g1; dst[5] = g1;
        dst[6] = g2; dst[7] = g2; dst[8] = g2;
        dst[9] = g3; dst[10] = g3; dst[11] = g3;
    }
#elif defined(HAS_SIMD_NEON)
    // Process 8 pixels at a time using NEON
    for (; i + 7 < count; i += 8) {
        uint8x8_t g = vld1_u8(&gray[i]);
        uint8x8x3_t rgb_vec = {{g, g, g}};
        vst3_u8(&rgb[i * 3], rgb_vec);
    }
#endif
    // Scalar remainder
    for (; i < count; i++) {
        uint8_t val = gray[i];
        rgb[i * 3] = val;
        rgb[i * 3 + 1] = val;
        rgb[i * 3 + 2] = val;
    }
}

// SIMD-optimized stretch application for x86_64 (SSE2)
#ifdef HAS_SIMD_SSE2
static inline void apply_stretch_simd_sse2(
    const float* channel_data,
    uint8_t* output,
    size_t count,
    float native_shadows,
    float native_highlights,
    float k1,
    float k2,
    float midtones,
    int stride  // 1 for grayscale-to-mono, 3 for RGB
) {
    const __m128 v_shadows = _mm_set1_ps(native_shadows);
    const __m128 v_highlights = _mm_set1_ps(native_highlights);
    const __m128 v_k1 = _mm_set1_ps(k1);
    const __m128 v_k2 = _mm_set1_ps(k2);
    const __m128 v_midtones = _mm_set1_ps(midtones);
    const __m128 v_zero = _mm_setzero_ps();
    const __m128 v_255 = _mm_set1_ps(255.0f);

    size_t i = 0;
    // Process 4 pixels at a time
    for (; i + 3 < count; i += 4) {
        __m128 input = _mm_loadu_ps(&channel_data[i]);

        // Create masks for the conditions
        __m128 mask_low = _mm_cmplt_ps(input, v_shadows);
        __m128 mask_high = _mm_cmpge_ps(input, v_highlights);

        // Compute stretch for middle range
        __m128 input_floored = _mm_sub_ps(input, v_shadows);
        __m128 numerator = _mm_mul_ps(input_floored, v_k1);
        __m128 denominator = _mm_sub_ps(_mm_mul_ps(input_floored, v_k2), v_midtones);
        __m128 output_val = _mm_div_ps(numerator, denominator);

        // Apply conditions: 0 if < shadows, 255 if >= highlights, computed otherwise
        output_val = _mm_and_ps(output_val, _mm_andnot_ps(mask_low, _mm_andnot_ps(mask_high, _mm_castsi128_ps(_mm_set1_epi32(-1)))));
        output_val = _mm_or_ps(output_val, _mm_and_ps(v_255, mask_high));

        // Clamp to [0, 255]
        output_val = _mm_min_ps(_mm_max_ps(output_val, v_zero), v_255);

        // Convert to int and store
        __m128i output_int = _mm_cvtps_epi32(output_val);
        output_int = _mm_packs_epi32(output_int, output_int);
        output_int = _mm_packus_epi16(output_int, output_int);

        // Extract and write
        uint32_t packed = _mm_cvtsi128_si32(output_int);
        for (int j = 0; j < 4; j++) {
            output[(i + j) * stride] = (uint8_t)((packed >> (j * 8)) & 0xFF);
        }
    }

    // Handle remaining pixels
    for (; i < count; i++) {
        float input = channel_data[i];
        float out;

        if (input < native_shadows) {
            out = 0.0f;
        } else if (input >= native_highlights) {
            out = 255.0f;
        } else {
            float input_floored = input - native_shadows;
            out = (input_floored * k1) / (input_floored * k2 - midtones);
        }

        output[i * stride] = (uint8_t)fmaxf(0.0f, fminf(255.0f, out));
    }
}
#endif

// SIMD-optimized stretch application for ARM (NEON)
#ifdef HAS_SIMD_NEON
static inline void apply_stretch_simd_neon(
    const float* channel_data,
    uint8_t* output,
    size_t count,
    float native_shadows,
    float native_highlights,
    float k1,
    float k2,
    float midtones,
    int stride
) {
    const float32x4_t v_shadows = vdupq_n_f32(native_shadows);
    const float32x4_t v_highlights = vdupq_n_f32(native_highlights);
    const float32x4_t v_k1 = vdupq_n_f32(k1);
    const float32x4_t v_k2 = vdupq_n_f32(k2);
    const float32x4_t v_midtones = vdupq_n_f32(midtones);
    const float32x4_t v_zero = vdupq_n_f32(0.0f);
    const float32x4_t v_255 = vdupq_n_f32(255.0f);

    size_t i = 0;
    // Process 4 pixels at a time
    for (; i + 3 < count; i += 4) {
        float32x4_t input = vld1q_f32(&channel_data[i]);

        // Create masks for the conditions
        uint32x4_t mask_low = vcltq_f32(input, v_shadows);
        uint32x4_t mask_high = vcgeq_f32(input, v_highlights);

        // Compute stretch for middle range
        float32x4_t input_floored = vsubq_f32(input, v_shadows);
        float32x4_t numerator = vmulq_f32(input_floored, v_k1);
        float32x4_t denominator = vsubq_f32(vmulq_f32(input_floored, v_k2), v_midtones);

        // Division using reciprocal estimate + Newton-Raphson refinement for better accuracy
        float32x4_t recip = vrecpeq_f32(denominator);
        recip = vmulq_f32(vrecpsq_f32(denominator, recip), recip);  // One N-R iteration
        float32x4_t output_val = vmulq_f32(numerator, recip);

        // Apply conditions using NEON select
        output_val = vbslq_f32(mask_low, v_zero, output_val);
        output_val = vbslq_f32(mask_high, v_255, output_val);

        // Clamp to [0, 255]
        output_val = vminq_f32(vmaxq_f32(output_val, v_zero), v_255);

        // Convert to int
        uint32x4_t output_int = vcvtq_u32_f32(output_val);

        // Narrow to 16-bit then 8-bit
        uint16x4_t output_16 = vmovn_u32(output_int);
        uint8x8_t output_8 = vmovn_u16(vcombine_u16(output_16, output_16));

        // Extract and write
        uint32_t packed[2];
        vst1_u8((uint8_t*)packed, output_8);
        for (int j = 0; j < 4; j++) {
            output[(i + j) * stride] = ((uint8_t*)packed)[j];
        }
    }

    // Handle remaining pixels
    for (; i < count; i++) {
        float input = channel_data[i];
        float out;

        if (input < native_shadows) {
            out = 0.0f;
        } else if (input >= native_highlights) {
            out = 255.0f;
        } else {
            float input_floored = input - native_shadows;
            out = (input_floored * k1) / (input_floored * k2 - midtones);
        }

        output[i * stride] = (uint8_t)fmaxf(0.0f, fminf(255.0f, out));
    }
}
#endif

// Generic wrapper that chooses the best SIMD implementation
static inline void apply_stretch_simd(
    const float* channel_data,
    uint8_t* output,
    size_t count,
    float native_shadows,
    float native_highlights,
    float k1,
    float k2,
    float midtones,
    int stride
) {
#ifdef HAS_SIMD_SSE2
    apply_stretch_simd_sse2(channel_data, output, count, native_shadows, native_highlights, k1, k2, midtones, stride);
#elif defined(HAS_SIMD_NEON)
    apply_stretch_simd_neon(channel_data, output, count, native_shadows, native_highlights, k1, k2, midtones, stride);
#else
    // Fallback scalar implementation
    for (size_t i = 0; i < count; i++) {
        float input = channel_data[i];
        float out;

        if (input < native_shadows) {
            out = 0.0f;
        } else if (input >= native_highlights) {
            out = 255.0f;
        } else {
            float input_floored = input - native_shadows;
            out = (input_floored * k1) / (input_floored * k2 - midtones);
        }

        output[i * stride] = (uint8_t)fmaxf(0.0f, fminf(255.0f, out));
    }
#endif
}

int process_fits_file(const char* fits_path, const ProcessConfig* config, ProcessedImage* out_image) {
    FitsMetadata meta;
    uint16_t* u16_data = NULL;
    float* f32_data = NULL;

    memset(out_image, 0, sizeof(ProcessedImage));

    // Read FITS image using custom reader
    if (read_fits_image(fits_path, &meta, &u16_data, &f32_data) != 0) {
        return -1;
    }

    // For now, only handle mono images
    if (meta.channels != 1) {
        set_error("Only mono images supported for now");
        if (u16_data) free(u16_data);
        if (f32_data) free(f32_data);
        return -1;
    }

    size_t width = meta.width;
    size_t height = meta.height;
    float* float_data = NULL;

    if (u16_data) {
        // Apply downscaling if needed
        if (config->downscale_factor > 1) {
            size_t new_w = width / config->downscale_factor;
            size_t new_h = height / config->downscale_factor;
            uint16_t* downscaled = malloc(new_w * new_h * sizeof(uint16_t));
            downscale_u16(u16_data, downscaled, width, height, config->downscale_factor);
            free(u16_data);
            u16_data = downscaled;
            width = new_w;
            height = new_h;
        }

        // Check for Bayer pattern and debayer
        if (config->apply_debayer && meta.bayer_pattern != BAYER_NONE) {
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
    } else if (f32_data) {
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
    } else {
        set_error("No pixel data read from FITS");
        return -1;
    }

    // Apply stretch and convert to 8-bit
    int num_channels = out_image->is_color ? 3 : 1;
    out_image->width = width;
    out_image->height = height;
    out_image->data = malloc(width * height * 3);  // Always RGB output

    if (config->auto_stretch) {
        // Process each channel (sequential without OpenMP on this platform)
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
                // Write to separate RGB channels using SIMD
                apply_stretch_simd(channel_data, &out_image->data[c], channel_size,
                                  native_shadows, native_highlights, k1, k2, midtones, 3);
            } else {
                // Grayscale - use SIMD then replicate to RGB
                uint8_t* temp = malloc(channel_size);
                apply_stretch_simd(channel_data, temp, channel_size,
                                  native_shadows, native_highlights, k1, k2, midtones, 1);
                // Replicate to RGB using optimized function
                replicate_gray_to_rgb(temp, out_image->data, channel_size);
                free(temp);
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

// 2x2 binning (averaging) for mono images in preview mode (SIMD optimized)
static void bin_2x2_float(const float* in, float* out, size_t w, size_t h) {
    size_t out_w = w / 2;
    size_t out_h = h / 2;

    for (size_t y = 0; y < out_h; y++) {
        size_t in_y = y * 2;
        const float* row0 = &in[in_y * w];
        const float* row1 = &in[(in_y + 1) * w];
        float* out_row = &out[y * out_w];
        size_t x = 0;

#ifdef HAS_SIMD_SSE2
        // Process 4 output pixels at a time (8 input pixels per row)
        __m128 quarter = _mm_set1_ps(0.25f);
        for (; x + 3 < out_w; x += 4) {
            size_t in_x = x * 2;
            // Load 8 floats from each row
            __m128 r0_a = _mm_loadu_ps(&row0[in_x]);     // p00,p01,p02,p03
            __m128 r0_b = _mm_loadu_ps(&row0[in_x + 4]); // p04,p05,p06,p07
            __m128 r1_a = _mm_loadu_ps(&row1[in_x]);
            __m128 r1_b = _mm_loadu_ps(&row1[in_x + 4]);

            // Horizontal pairwise add: (a+b, c+d, e+f, g+h)
            __m128 sum0 = _mm_add_ps(_mm_shuffle_ps(r0_a, r0_b, _MM_SHUFFLE(2,0,2,0)),
                                     _mm_shuffle_ps(r0_a, r0_b, _MM_SHUFFLE(3,1,3,1)));
            __m128 sum1 = _mm_add_ps(_mm_shuffle_ps(r1_a, r1_b, _MM_SHUFFLE(2,0,2,0)),
                                     _mm_shuffle_ps(r1_a, r1_b, _MM_SHUFFLE(3,1,3,1)));

            // Add rows and multiply by 0.25
            __m128 result = _mm_mul_ps(_mm_add_ps(sum0, sum1), quarter);
            _mm_storeu_ps(&out_row[x], result);
        }
#elif defined(HAS_SIMD_NEON)
        float32x4_t quarter = vdupq_n_f32(0.25f);
        for (; x + 3 < out_w; x += 4) {
            size_t in_x = x * 2;
            float32x4x2_t r0 = vld2q_f32(&row0[in_x]);  // Deinterleave pairs
            float32x4x2_t r1 = vld2q_f32(&row1[in_x]);

            float32x4_t sum0 = vaddq_f32(r0.val[0], r0.val[1]);  // p00+p01, p02+p03, ...
            float32x4_t sum1 = vaddq_f32(r1.val[0], r1.val[1]);

            float32x4_t result = vmulq_f32(vaddq_f32(sum0, sum1), quarter);
            vst1q_f32(&out_row[x], result);
        }
#endif
        // Scalar remainder
        for (; x < out_w; x++) {
            size_t in_x = x * 2;
            float p00 = row0[in_x];
            float p01 = row0[in_x + 1];
            float p10 = row1[in_x];
            float p11 = row1[in_x + 1];
            out_row[x] = (p00 + p01 + p10 + p11) * 0.25f;
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

    // Compute median using quickselect (faster than full sort)
    float median = find_median(samples, num_samples);

    // Compute MADN (Median Absolute Deviation Normalized)
    float* deviations = malloc(num_samples * sizeof(float));
    for (size_t i = 0; i < num_samples; i++) {
        deviations[i] = fabsf(samples[i] - median);
    }
    float madn = 1.4826f * find_median(deviations, num_samples);

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

}

void free_processed_image(ProcessedImage* image) {
    if (image && image->data) {
        free(image->data);
        image->data = NULL;
    }
}

// Process XISF file using the same pipeline as FITS
static int process_xisf_internal(const char* xisf_path, const ProcessConfig* config, ProcessedImage* out_image) {
    FitsMetadata meta;
    uint16_t* u16_data = NULL;
    float* float_data = NULL;

    memset(out_image, 0, sizeof(ProcessedImage));

    // Read XISF image (returns float32 data)
    if (read_xisf_image(xisf_path, &meta, &u16_data, &float_data) != 0) {
        return -1;
    }

    // XISF reader always returns float data
    if (!float_data) {
        set_error("XISF reader returned no data");
        return -1;
    }

    size_t width = meta.width;
    size_t height = meta.height;

    // Apply downscaling if needed
    if (config->downscale_factor > 1) {
        size_t new_w = width / config->downscale_factor;
        size_t new_h = height / config->downscale_factor;
        size_t num_channels = meta.channels;
        float* downscaled = malloc(new_w * new_h * num_channels * sizeof(float));

        // Downscale each channel
        for (size_t c = 0; c < num_channels; c++) {
            float* src_plane = float_data + c * width * height;
            float* dst_plane = downscaled + c * new_w * new_h;
            for (size_t y = 0; y < new_h; y++) {
                for (size_t x = 0; x < new_w; x++) {
                    dst_plane[y * new_w + x] = src_plane[(y * config->downscale_factor) * width + (x * config->downscale_factor)];
                }
            }
        }

        free(float_data);
        float_data = downscaled;
        width = new_w;
        height = new_h;
    }

    // Handle mono images with Bayer pattern
    if (meta.channels == 1 && config->apply_debayer && meta.bayer_pattern != BAYER_NONE) {
        size_t out_w = width / 2;
        size_t out_h = height / 2;
        float* rgb_data = malloc(out_w * out_h * 3 * sizeof(float));
        super_pixel_debayer_f32(float_data, rgb_data, width, height, meta.bayer_pattern);
        free(float_data);
        float_data = rgb_data;
        width = out_w;
        height = out_h;
        out_image->is_color = 1;
        meta.channels = 3;
    } else if (meta.channels == 3) {
        // Already RGB
        out_image->is_color = 1;
    } else {
        // Apply 2x2 binning in preview mode for mono images
        if (config->preview_mode && meta.channels == 1 && meta.bayer_pattern == BAYER_NONE) {
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

            float max_input = 65536.0f;  // Standard 16-bit range

            compute_stretch_params(channel_data, channel_size, max_input,
                                 &shadows, &highlights, &midtones);

            float hs_range_factor = (highlights == shadows) ? 1.0f : 1.0f / (highlights - shadows);
            float native_shadows = shadows * max_input;
            float native_highlights = highlights * max_input;
            float k1 = (midtones - 1.0f) * hs_range_factor * 255.0f / max_input;
            float k2 = ((2.0f * midtones) - 1.0f) * hs_range_factor / max_input;

            if (out_image->is_color) {
                apply_stretch_simd(channel_data, &out_image->data[c], channel_size,
                                  native_shadows, native_highlights, k1, k2, midtones, 3);
            } else {
                uint8_t* temp = malloc(channel_size);
                apply_stretch_simd(channel_data, temp, channel_size,
                                  native_shadows, native_highlights, k1, k2, midtones, 1);
                replicate_gray_to_rgb(temp, out_image->data, channel_size);
                free(temp);
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

// Main entry point - handles both FITS and XISF files
int process_image_file(const char* path, const ProcessConfig* config, ProcessedImage* out_image) {
    if (is_xisf_file(path)) {
        return process_xisf_internal(path, config, out_image);
    } else {
        return process_fits_file(path, config, out_image);
    }
}
