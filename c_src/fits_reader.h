#ifndef FITS_READER_H
#define FITS_READER_H

#include "fits_processor.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// FITS constants
#define FITS_BLOCK_SIZE 2880
#define FITS_CARD_SIZE 80

/**
 * Check if a file is a FITS file (by extension or magic bytes)
 * @param path File path
 * @return 1 if FITS, 0 otherwise
 */
int is_fits_file(const char* path);

/**
 * Read FITS image and return pixel data compatible with processing pipeline
 *
 * @param fits_path Path to FITS file
 * @param meta Output metadata (FitsMetadata)
 * @param out_data_u16 Output uint16 data (if BITPIX=16)
 * @param out_data_f32 Output float32 data (if BITPIX=-32)
 * @return 0 on success, -1 on error
 */
int read_fits_image(
    const char* fits_path,
    FitsMetadata* meta,
    uint16_t** out_data_u16,
    float** out_data_f32
);

#ifdef __cplusplus
}
#endif

#endif // FITS_READER_H
