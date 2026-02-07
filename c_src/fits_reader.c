#include "fits_reader.h"
#include "compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

// Thread-local error buffer (shared with fits_processor.c)
extern THREAD_LOCAL char error_buffer[512];

static void set_fits_error(const char* msg) {
    snprintf(error_buffer, sizeof(error_buffer), "%s", msg);
}

// FITS uses big-endian byte order - use compiler builtins for single-instruction swaps
static inline uint16_t swap16(uint16_t val) {
    return BSWAP16(val);
}

static inline int16_t swap16s(int16_t val) {
    return (int16_t)BSWAP16((uint16_t)val);
}

static inline uint32_t swap32(uint32_t val) {
    return BSWAP32(val);
}

static inline uint64_t swap64(uint64_t val) {
    return BSWAP64(val);
}

static inline float swap_float(float val) {
    uint32_t tmp;
    memcpy(&tmp, &val, 4);
    tmp = BSWAP32(tmp);
    float result;
    memcpy(&result, &tmp, 4);
    return result;
}

int is_fits_file(const char* path) {
    // Check extension
    const char* ext = strrchr(path, '.');
    if (ext) {
        if (STRCASECMP(ext, ".fits") == 0 || STRCASECMP(ext, ".fit") == 0) {
            return 1;
        }
    }

    // Check magic bytes: "SIMPLE  ="
    FILE* f = fopen(path, "rb");
    if (!f) return 0;

    char buf[10];
    size_t read = fread(buf, 1, 9, f);
    fclose(f);

    if (read >= 9 && strncmp(buf, "SIMPLE  =", 9) == 0) {
        return 1;
    }

    return 0;
}

// Parse a single keyword value from an 80-char card
// Returns pointer to value start, or NULL if not found
static const char* get_keyword_value(const char* card, const char* keyword) {
    size_t kw_len = strlen(keyword);
    if (strncmp(card, keyword, kw_len) != 0) return NULL;

    // Find '=' sign
    const char* eq = strchr(card, '=');
    if (!eq) return NULL;

    // Skip whitespace after '='
    const char* val = eq + 1;
    while (*val == ' ' && val < card + 80) val++;

    return val;
}

// Parse integer value from card
static int parse_int_keyword(const char* card, const char* keyword, int* value) {
    const char* val = get_keyword_value(card, keyword);
    if (!val) return -1;
    *value = atoi(val);
    return 0;
}

// Parse float value from card
static int parse_float_keyword(const char* card, const char* keyword, double* value) {
    const char* val = get_keyword_value(card, keyword);
    if (!val) return -1;
    *value = atof(val);
    return 0;
}

// Parse string value from card (handles 'quoted strings')
static int parse_string_keyword(const char* card, const char* keyword, char* value, size_t max_len) {
    const char* val = get_keyword_value(card, keyword);
    if (!val) return -1;

    // Skip leading quote
    if (*val == '\'') val++;

    // Copy until closing quote or end
    size_t i = 0;
    while (*val && *val != '\'' && i < max_len - 1 && val < card + 80) {
        value[i++] = *val++;
    }
    value[i] = '\0';

    // Trim trailing spaces
    while (i > 0 && value[i-1] == ' ') {
        value[--i] = '\0';
    }

    return 0;
}

// FITS header info
typedef struct {
    int bitpix;
    int naxis;
    int naxis1;  // width
    int naxis2;  // height
    int naxis3;  // channels (for 3D)
    double bzero;
    double bscale;
    char bayerpat[32];
    char roworder[32];
    size_t header_size;  // Total header size in bytes
} FitsHeader;

// Read and parse FITS header
static int read_fits_header(FILE* f, FitsHeader* hdr) {
    memset(hdr, 0, sizeof(FitsHeader));
    hdr->bscale = 1.0;  // Default
    hdr->bzero = 0.0;   // Default

    char block[FITS_BLOCK_SIZE];
    int found_end = 0;
    size_t total_read = 0;

    while (!found_end) {
        size_t read = fread(block, 1, FITS_BLOCK_SIZE, f);
        if (read != FITS_BLOCK_SIZE) {
            set_fits_error("Failed to read FITS header block");
            return -1;
        }
        total_read += read;

        // Parse each 80-char card in this block
        for (int i = 0; i < FITS_BLOCK_SIZE / FITS_CARD_SIZE; i++) {
            char* card = block + i * FITS_CARD_SIZE;

            // Check for END keyword
            if (strncmp(card, "END", 3) == 0 && (card[3] == ' ' || card[3] == '\0')) {
                found_end = 1;
                break;
            }

            // Parse keywords we care about
            int ival;
            double dval;
            char sval[80];

            if (parse_int_keyword(card, "BITPIX  ", &ival) == 0) {
                hdr->bitpix = ival;
            } else if (parse_int_keyword(card, "NAXIS   ", &ival) == 0) {
                hdr->naxis = ival;
            } else if (parse_int_keyword(card, "NAXIS1", &ival) == 0) {
                hdr->naxis1 = ival;
            } else if (parse_int_keyword(card, "NAXIS2", &ival) == 0) {
                hdr->naxis2 = ival;
            } else if (parse_int_keyword(card, "NAXIS3", &ival) == 0) {
                hdr->naxis3 = ival;
            } else if (parse_float_keyword(card, "BZERO", &dval) == 0) {
                hdr->bzero = dval;
            } else if (parse_float_keyword(card, "BSCALE", &dval) == 0) {
                hdr->bscale = dval;
            } else if (parse_string_keyword(card, "BAYERPAT", sval, sizeof(sval)) == 0) {
                strncpy(hdr->bayerpat, sval, sizeof(hdr->bayerpat) - 1);
            } else if (parse_string_keyword(card, "ROWORDER", sval, sizeof(sval)) == 0) {
                strncpy(hdr->roworder, sval, sizeof(hdr->roworder) - 1);
            }
        }
    }

    hdr->header_size = total_read;

    // Validate required fields
    if (hdr->bitpix == 0) {
        set_fits_error("Missing BITPIX keyword in FITS header");
        return -1;
    }
    if (hdr->naxis < 2) {
        set_fits_error("FITS image must have at least 2 dimensions");
        return -1;
    }
    if (hdr->naxis1 <= 0 || hdr->naxis2 <= 0) {
        set_fits_error("Invalid FITS image dimensions");
        return -1;
    }

    return 0;
}

int read_fits_image(const char* fits_path, FitsMetadata* meta,
                    uint16_t** out_data_u16, float** out_data_f32) {
    FILE* f = fopen(fits_path, "rb");
    if (!f) {
        set_fits_error("Failed to open FITS file");
        return -1;
    }

    // Parse header
    FitsHeader hdr;
    if (read_fits_header(f, &hdr) != 0) {
        fclose(f);
        return -1;
    }

    // Fill metadata
    meta->width = hdr.naxis1;
    meta->height = hdr.naxis2;
    meta->channels = (hdr.naxis >= 3 && hdr.naxis3 > 0) ? hdr.naxis3 : 1;

    // Parse Bayer pattern
    meta->bayer_pattern = BAYER_NONE;
    if (hdr.bayerpat[0]) {
        if (strcmp(hdr.bayerpat, "RGGB") == 0) meta->bayer_pattern = BAYER_RGGB;
        else if (strcmp(hdr.bayerpat, "BGGR") == 0) meta->bayer_pattern = BAYER_BGGR;
        else if (strcmp(hdr.bayerpat, "GBRG") == 0) meta->bayer_pattern = BAYER_GBRG;
        else if (strcmp(hdr.bayerpat, "GRBG") == 0) meta->bayer_pattern = BAYER_GRBG;
    }

    // Parse row order
    meta->flip_vertical = (strcmp(hdr.roworder, "TOP-DOWN") == 0) ? 1 : 0;

    // Calculate data size
    size_t num_pixels = meta->width * meta->height * meta->channels;
    int bytes_per_pixel = abs(hdr.bitpix) / 8;
    size_t data_size = num_pixels * bytes_per_pixel;

    // Read raw data
    uint8_t* raw_data = malloc(data_size);
    if (!raw_data) {
        set_fits_error("Out of memory for FITS data");
        fclose(f);
        return -1;
    }

    size_t read = fread(raw_data, 1, data_size, f);
    fclose(f);

    if (read != data_size) {
        char err[256];
        snprintf(err, sizeof(err), "Failed to read FITS data: expected %zu, got %zu", data_size, read);
        set_fits_error(err);
        free(raw_data);
        return -1;
    }

    // Convert based on BITPIX
    *out_data_u16 = NULL;
    *out_data_f32 = NULL;

    if (hdr.bitpix == 16) {
        // 16-bit signed integer (often with BZERO=32768 for unsigned)
        meta->dtype = DTYPE_UINT16;
        uint16_t* u16_data = malloc(num_pixels * sizeof(uint16_t));
        if (!u16_data) {
            set_fits_error("Out of memory");
            free(raw_data);
            return -1;
        }

        uint16_t* src = (uint16_t*)raw_data;

        // Fast path: BZERO=32768, BSCALE=1 (common for unsigned 16-bit)
        // This converts signed to unsigned by adding 32768, which is just XOR with 0x8000
        if (hdr.bzero == 32768.0 && hdr.bscale == 1.0) {
            for (size_t i = 0; i < num_pixels; i++) {
                u16_data[i] = swap16(src[i]) ^ 0x8000;
            }
        } else if (hdr.bzero == 0.0 && hdr.bscale == 1.0) {
            // Simple swap only - no scaling needed
            for (size_t i = 0; i < num_pixels; i++) {
                u16_data[i] = swap16(src[i]);
            }
        } else {
            // General case with BZERO/BSCALE
            int16_t* src_s = (int16_t*)raw_data;
            for (size_t i = 0; i < num_pixels; i++) {
                int16_t val = swap16s(src_s[i]);
                double scaled = hdr.bzero + hdr.bscale * val;
                if (scaled < 0) scaled = 0;
                if (scaled > 65535) scaled = 65535;
                u16_data[i] = (uint16_t)scaled;
            }
        }

        free(raw_data);
        *out_data_u16 = u16_data;

    } else if (hdr.bitpix == -32) {
        // 32-bit float
        meta->dtype = DTYPE_FLOAT32;
        float* f32_data = malloc(num_pixels * sizeof(float));
        if (!f32_data) {
            set_fits_error("Out of memory");
            free(raw_data);
            return -1;
        }

        float* src = (float*)raw_data;
        for (size_t i = 0; i < num_pixels; i++) {
            // Swap bytes and apply BZERO/BSCALE
            float val = swap_float(src[i]);
            f32_data[i] = (float)(hdr.bzero + hdr.bscale * val);
        }

        free(raw_data);
        *out_data_f32 = f32_data;

    } else if (hdr.bitpix == 8) {
        // 8-bit unsigned - convert to uint16
        meta->dtype = DTYPE_UINT16;
        uint16_t* u16_data = malloc(num_pixels * sizeof(uint16_t));
        if (!u16_data) {
            set_fits_error("Out of memory");
            free(raw_data);
            return -1;
        }

        for (size_t i = 0; i < num_pixels; i++) {
            double scaled = hdr.bzero + hdr.bscale * raw_data[i];
            u16_data[i] = (uint16_t)(scaled * 256);  // Scale 8-bit to 16-bit range
        }

        free(raw_data);
        *out_data_u16 = u16_data;

    } else if (hdr.bitpix == 32) {
        // 32-bit integer - convert to float
        meta->dtype = DTYPE_FLOAT32;
        float* f32_data = malloc(num_pixels * sizeof(float));
        if (!f32_data) {
            set_fits_error("Out of memory");
            free(raw_data);
            return -1;
        }

        int32_t* src = (int32_t*)raw_data;
        for (size_t i = 0; i < num_pixels; i++) {
            int32_t val = (int32_t)swap32((uint32_t)src[i]);
            f32_data[i] = (float)(hdr.bzero + hdr.bscale * val);
        }

        free(raw_data);
        *out_data_f32 = f32_data;

    } else if (hdr.bitpix == -64) {
        // 64-bit float - convert to float32
        meta->dtype = DTYPE_FLOAT32;
        float* f32_data = malloc(num_pixels * sizeof(float));
        if (!f32_data) {
            set_fits_error("Out of memory");
            free(raw_data);
            return -1;
        }

        uint64_t* src = (uint64_t*)raw_data;
        for (size_t i = 0; i < num_pixels; i++) {
            uint64_t val = swap64(src[i]);
            double dval;
            memcpy(&dval, &val, 8);
            f32_data[i] = (float)(hdr.bzero + hdr.bscale * dval);
        }

        free(raw_data);
        *out_data_f32 = f32_data;

    } else {
        char err[256];
        snprintf(err, sizeof(err), "Unsupported BITPIX value: %d", hdr.bitpix);
        set_fits_error(err);
        free(raw_data);
        return -1;
    }

    return 0;
}
