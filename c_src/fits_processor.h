#ifndef FITS_PROCESSOR_H
#define FITS_PROCESSOR_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Bayer pattern types
typedef enum {
    BAYER_NONE = 0,
    BAYER_RGGB = 1,
    BAYER_BGGR = 2,
    BAYER_GBRG = 3,
    BAYER_GRBG = 4
} BayerPattern;

// Image data type
typedef enum {
    DTYPE_UINT16 = 0,
    DTYPE_FLOAT32 = 1
} DataType;

// FITS image metadata
typedef struct {
    size_t width;
    size_t height;
    size_t channels;
    DataType dtype;
    BayerPattern bayer_pattern;
    int flip_vertical;  // 1 if ROWORDER is TOP-DOWN
} FitsMetadata;

// Processed image output (always 8-bit RGB)
typedef struct {
    uint8_t* data;      // RGB data (width * height * 3)
    size_t width;
    size_t height;
    int is_color;       // 1 if RGB, 0 if grayscale
} ProcessedImage;

// Processing configuration
typedef struct {
    int downscale_factor;   // 1 = no downscaling
    int jpeg_quality;       // 1-100
    int apply_debayer;      // 1 = apply debayering if Bayer pattern detected
    int preview_mode;       // 1 = apply 2x2 binning for mono images (fast preview)
    int auto_stretch;       // 1 = auto-stretch, 0 = use manual params
    float manual_shadows;   // Manual stretch parameters (if auto_stretch=0)
    float manual_highlights;
    float manual_midtones;
} ProcessConfig;

/**
 * Process FITS file and return RGB image data
 *
 * @param fits_path Path to FITS file
 * @param config Processing configuration
 * @param out_image Output image (caller must free with free_processed_image)
 * @return 0 on success, -1 on error
 */
int process_fits_file(
    const char* fits_path,
    const ProcessConfig* config,
    ProcessedImage* out_image
);

/**
 * Save processed image as JPEG
 *
 * @param image Processed image
 * @param output_path Output JPEG path
 * @param quality JPEG quality (1-100)
 * @return 0 on success, -1 on error
 */
int save_jpeg(
    const ProcessedImage* image,
    const char* output_path,
    int quality
);

/**
 * Free processed image memory
 */
void free_processed_image(ProcessedImage* image);

/**
 * Get last error message
 */
const char* get_last_error(void);

#ifdef __cplusplus
}
#endif

#endif // FITS_PROCESSOR_H
