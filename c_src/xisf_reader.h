#ifndef XISF_READER_H
#define XISF_READER_H

#include "fits_processor.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// XISF block location types
typedef enum {
    XISF_LOCATION_ATTACHMENT = 0,  // Raw bytes at file offset
    XISF_LOCATION_EMBEDDED = 1,    // Base64 encoded in XML
    XISF_LOCATION_INLINE = 2       // Hex encoded in XML
} XisfBlockLocation;

// XISF compression types
typedef enum {
    XISF_COMPRESS_NONE = 0,
    XISF_COMPRESS_ZLIB = 1,
    XISF_COMPRESS_LZ4 = 2,
    XISF_COMPRESS_LZ4HC = 3,
    XISF_COMPRESS_ZSTD = 4
} XisfCompression;

// XISF sample formats
typedef enum {
    XISF_SAMPLE_UINT8 = 0,
    XISF_SAMPLE_UINT16 = 1,
    XISF_SAMPLE_UINT32 = 2,
    XISF_SAMPLE_FLOAT32 = 3,
    XISF_SAMPLE_FLOAT64 = 4
} XisfSampleFormat;

// XISF image metadata parsed from XML
typedef struct {
    size_t width;
    size_t height;
    size_t channels;
    XisfSampleFormat sample_format;
    int is_planar;           // 1 if planar storage (RRRGGGBBB), 0 if normal (RGBRGBRGB)

    // Block location info
    XisfBlockLocation location;
    uint64_t attachment_pos;   // File offset for attachment
    uint64_t block_size;       // Size of (possibly compressed) data

    // Compression info
    XisfCompression compression;
    uint64_t uncompressed_size;
    int byte_shuffled;      // 1 if data was byte-shuffled before compression
    int shuffle_item_size;  // Item size for byte shuffling (e.g., 2 for UInt16)

    // Optional: embedded/inline data pointer (points into XML buffer)
    const char* embedded_data;
    size_t embedded_len;

    // Bayer pattern (if present in XISF properties)
    BayerPattern bayer_pattern;

    // Row order
    int flip_vertical;
} XisfImageInfo;

/**
 * Check if a file is an XISF file (by extension or magic bytes)
 * @param path File path
 * @return 1 if XISF, 0 otherwise
 */
int is_xisf_file(const char* path);

/**
 * Process XISF file and return pixel data compatible with FITS pipeline
 *
 * @param xisf_path Path to XISF file
 * @param meta Output metadata (FitsMetadata for pipeline compatibility)
 * @param out_data_u16 Output uint16 data (if sample format is UInt8/UInt16)
 * @param out_data_f32 Output float32 data (if sample format is Float32/Float64)
 * @return 0 on success, -1 on error
 */
int read_xisf_image(
    const char* xisf_path,
    FitsMetadata* meta,
    uint16_t** out_data_u16,
    float** out_data_f32
);

/**
 * Base64 decode helper
 * @param input Base64 encoded string
 * @param input_len Length of input
 * @param output Output buffer (caller allocates)
 * @param output_len Actual decoded length (output)
 * @return 0 on success, -1 on error
 */
int base64_decode(const char* input, size_t input_len, uint8_t* output, size_t* output_len);

/**
 * Get decoded size from base64 string length
 */
size_t base64_decoded_size(size_t input_len);

#ifdef __cplusplus
}
#endif

#endif // XISF_READER_H
