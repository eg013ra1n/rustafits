#ifndef FITS_CONVERTER_H
#define FITS_CONVERTER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Opaque handle to FitsConverter
typedef struct FitsConverterHandle FitsConverterHandle;

// Error codes
#define FITS_OK 0
#define FITS_ERROR_NULL_POINTER -1
#define FITS_ERROR_INVALID_UTF8 -2
#define FITS_ERROR_CONVERSION -3

/**
 * Create a new FITS converter with default settings
 * @return Opaque pointer to FitsConverter (NULL on error)
 */
FitsConverterHandle* fits_converter_create(void);

/**
 * Set downscale factor
 * @param handle Converter handle
 * @param factor Downscale factor (1 = no downscaling)
 * @return FITS_OK on success
 */
int fits_converter_set_downscale(FitsConverterHandle* handle, size_t factor);

/**
 * Set JPEG quality
 * @param handle Converter handle
 * @param quality JPEG quality (1-100)
 * @return FITS_OK on success
 */
int fits_converter_set_quality(FitsConverterHandle* handle, int quality);

/**
 * Disable automatic debayering
 * @param handle Converter handle
 * @return FITS_OK on success
 */
int fits_converter_disable_debayer(FitsConverterHandle* handle);

/**
 * Convert FITS file to JPEG
 * @param handle Converter handle
 * @param input_path Path to input FITS file (null-terminated C string)
 * @param output_path Path to output JPEG file (null-terminated C string)
 * @return FITS_OK on success, error code on failure
 */
int fits_converter_convert(
    FitsConverterHandle* handle,
    const char* input_path,
    const char* output_path
);

/**
 * Get last error message
 * @return Pointer to error message string (may be NULL)
 */
const char* fits_converter_get_error(void);

/**
 * Destroy the FITS converter and free memory
 * @param handle Converter handle to destroy
 */
void fits_converter_destroy(FitsConverterHandle* handle);

#ifdef __cplusplus
}
#endif

#endif // FITS_CONVERTER_H
