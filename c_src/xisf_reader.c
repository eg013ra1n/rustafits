#include "xisf_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <lz4.h>
#include <zstd.h>

// SIMD support detection
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #include <emmintrin.h>  // SSE2
    #define HAS_SIMD_SSE2 1
#elif defined(__ARM_NEON) || defined(__aarch64__)
    #include <arm_neon.h>
    #define HAS_SIMD_NEON 1
#endif

// Thread-local error buffer (shared with fits_processor.c via extern)
extern __thread char error_buffer[512];

static void set_xisf_error(const char* msg) {
    snprintf(error_buffer, sizeof(error_buffer), "%s", msg);
}

// XISF magic signature: "XISF0100"
static const char XISF_SIGNATURE[8] = {'X', 'I', 'S', 'F', '0', '1', '0', '0'};

int is_xisf_file(const char* path) {
    // Check extension first
    const char* ext = strrchr(path, '.');
    if (ext && strcasecmp(ext, ".xisf") == 0) {
        return 1;
    }

    // Check magic bytes
    FILE* f = fopen(path, "rb");
    if (!f) return 0;

    char sig[8];
    size_t read = fread(sig, 1, 8, f);
    fclose(f);

    if (read == 8 && memcmp(sig, XISF_SIGNATURE, 8) == 0) {
        return 1;
    }

    return 0;
}

// Helper: find attribute value in XML string
// Returns pointer to start of value (after opening quote), sets end to closing quote
static const char* find_xml_attr(const char* xml, const char* attr_name, const char** end) {
    char pattern[128];
    snprintf(pattern, sizeof(pattern), "%s=\"", attr_name);

    const char* pos = strstr(xml, pattern);
    if (!pos) {
        snprintf(pattern, sizeof(pattern), "%s='", attr_name);
        pos = strstr(xml, pattern);
    }
    if (!pos) return NULL;

    pos += strlen(pattern);
    char quote = *(pos - 1);  // " or '

    const char* close = strchr(pos, quote);
    if (!close) return NULL;

    *end = close;
    return pos;
}

// Helper: extract string attribute
static int get_xml_attr_str(const char* xml, const char* attr_name, char* out, size_t out_size) {
    const char* end;
    const char* start = find_xml_attr(xml, attr_name, &end);
    if (!start) return -1;

    size_t len = end - start;
    if (len >= out_size) len = out_size - 1;
    memcpy(out, start, len);
    out[len] = '\0';
    return 0;
}

// Parse geometry attribute "width:height:channels"
static int parse_geometry(const char* geom, size_t* width, size_t* height, size_t* channels) {
    char buf[128];
    strncpy(buf, geom, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';

    char* tok = strtok(buf, ":");
    if (!tok) return -1;
    *width = (size_t)atoll(tok);

    tok = strtok(NULL, ":");
    if (!tok) return -1;
    *height = (size_t)atoll(tok);

    tok = strtok(NULL, ":");
    if (tok) {
        *channels = (size_t)atoll(tok);
    } else {
        *channels = 1;  // Default to mono
    }

    return 0;
}

// Parse sample format string
static XisfSampleFormat parse_sample_format(const char* fmt) {
    if (strcasecmp(fmt, "UInt8") == 0) return XISF_SAMPLE_UINT8;
    if (strcasecmp(fmt, "UInt16") == 0) return XISF_SAMPLE_UINT16;
    if (strcasecmp(fmt, "UInt32") == 0) return XISF_SAMPLE_UINT32;
    if (strcasecmp(fmt, "Float32") == 0) return XISF_SAMPLE_FLOAT32;
    if (strcasecmp(fmt, "Float64") == 0) return XISF_SAMPLE_FLOAT64;
    return XISF_SAMPLE_FLOAT32;  // Default
}

// Parse location attribute "attachment:position:size" or "embedded" or "inline"
static int parse_location(const char* loc, XisfImageInfo* info) {
    if (strncmp(loc, "attachment:", 11) == 0) {
        info->location = XISF_LOCATION_ATTACHMENT;
        // Parse position:size
        char buf[128];
        strncpy(buf, loc + 11, sizeof(buf) - 1);
        buf[sizeof(buf) - 1] = '\0';

        char* colon = strchr(buf, ':');
        if (!colon) return -1;
        *colon = '\0';

        info->attachment_pos = (uint64_t)atoll(buf);
        info->block_size = (uint64_t)atoll(colon + 1);
        return 0;
    } else if (strcmp(loc, "embedded") == 0) {
        info->location = XISF_LOCATION_EMBEDDED;
        return 0;
    } else if (strcmp(loc, "inline") == 0) {
        info->location = XISF_LOCATION_INLINE;
        return 0;
    }
    return -1;
}

// Parse compression attribute "codec:uncompressed_size" or "codec:uncompressed_size:item_size"
// Codec can have +sh suffix for byte shuffling (e.g., "lz4hc+sh:51883264:2")
static int parse_compression(const char* comp, XisfImageInfo* info) {
    char buf[128];
    strncpy(buf, comp, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';

    char* colon = strchr(buf, ':');
    if (!colon) return -1;
    *colon = '\0';

    // Check for byte shuffling suffix (+sh)
    char* shuffle = strstr(buf, "+sh");
    if (shuffle) {
        *shuffle = '\0';  // Remove +sh from codec name
        info->byte_shuffled = 1;
    } else {
        info->byte_shuffled = 0;
    }

    // Parse codec
    if (strcasecmp(buf, "zlib") == 0) {
        info->compression = XISF_COMPRESS_ZLIB;
    } else if (strcasecmp(buf, "lz4") == 0) {
        info->compression = XISF_COMPRESS_LZ4;
    } else if (strcasecmp(buf, "lz4hc") == 0 || strcasecmp(buf, "lz4+hc") == 0) {
        info->compression = XISF_COMPRESS_LZ4HC;
    } else if (strcasecmp(buf, "zstd") == 0) {
        info->compression = XISF_COMPRESS_ZSTD;
    } else {
        return -1;  // Unknown codec
    }

    // Parse uncompressed size
    char* rest = colon + 1;
    char* next_colon = strchr(rest, ':');
    if (next_colon) {
        *next_colon = '\0';
        info->uncompressed_size = (uint64_t)atoll(rest);
        info->shuffle_item_size = (int)atoi(next_colon + 1);  // Item size for shuffling
    } else {
        info->uncompressed_size = (uint64_t)atoll(rest);
        info->shuffle_item_size = 1;
    }

    return 0;
}

// Parse XISF XML header to extract image info
static int parse_xisf_xml(const char* xml, size_t xml_len, XisfImageInfo* info) {
    memset(info, 0, sizeof(XisfImageInfo));
    info->bayer_pattern = BAYER_NONE;

    // Find <Image element
    const char* img_start = strstr(xml, "<Image");
    if (!img_start) {
        set_xisf_error("No <Image> element found in XISF header");
        return -1;
    }

    // Find end of Image element (either /> or >)
    const char* img_end = strchr(img_start, '>');
    if (!img_end) {
        set_xisf_error("Malformed <Image> element");
        return -1;
    }

    // Create a working copy of just the Image element
    size_t img_len = img_end - img_start + 1;
    char* img_elem = malloc(img_len + 1);
    memcpy(img_elem, img_start, img_len);
    img_elem[img_len] = '\0';

    // Parse geometry (required)
    char geom[128];
    if (get_xml_attr_str(img_elem, "geometry", geom, sizeof(geom)) != 0) {
        set_xisf_error("Missing geometry attribute in XISF Image");
        free(img_elem);
        return -1;
    }
    if (parse_geometry(geom, &info->width, &info->height, &info->channels) != 0) {
        set_xisf_error("Invalid geometry format");
        free(img_elem);
        return -1;
    }

    // Parse sampleFormat (default: Float32)
    char sample_fmt[32] = "Float32";
    get_xml_attr_str(img_elem, "sampleFormat", sample_fmt, sizeof(sample_fmt));
    info->sample_format = parse_sample_format(sample_fmt);

    // Parse colorSpace (for info, we handle Gray and RGB)
    char color_space[32] = "Gray";
    get_xml_attr_str(img_elem, "colorSpace", color_space, sizeof(color_space));
    // channels already parsed from geometry

    // Parse pixelStorage (normal or planar)
    char storage[32] = "planar";  // PixInsight default is planar
    get_xml_attr_str(img_elem, "pixelStorage", storage, sizeof(storage));
    info->is_planar = (strcasecmp(storage, "planar") == 0) ? 1 : 0;

    // Parse location (required)
    char location[256];
    if (get_xml_attr_str(img_elem, "location", location, sizeof(location)) != 0) {
        set_xisf_error("Missing location attribute in XISF Image");
        free(img_elem);
        return -1;
    }
    if (parse_location(location, info) != 0) {
        set_xisf_error("Invalid location format");
        free(img_elem);
        return -1;
    }

    // Parse compression (optional)
    char compression[128];
    if (get_xml_attr_str(img_elem, "compression", compression, sizeof(compression)) == 0) {
        if (parse_compression(compression, info) != 0) {
            char err[256];
            snprintf(err, sizeof(err), "Unsupported compression: %s", compression);
            set_xisf_error(err);
            free(img_elem);
            return -1;
        }
    } else {
        info->compression = XISF_COMPRESS_NONE;
        // For uncompressed, calculate expected size
        size_t bytes_per_sample = 0;
        switch (info->sample_format) {
            case XISF_SAMPLE_UINT8: bytes_per_sample = 1; break;
            case XISF_SAMPLE_UINT16: bytes_per_sample = 2; break;
            case XISF_SAMPLE_UINT32: bytes_per_sample = 4; break;
            case XISF_SAMPLE_FLOAT32: bytes_per_sample = 4; break;
            case XISF_SAMPLE_FLOAT64: bytes_per_sample = 8; break;
        }
        info->uncompressed_size = info->width * info->height * info->channels * bytes_per_sample;
    }

    free(img_elem);

    // Look for Bayer pattern in XISF properties
    const char* bayer_prop = strstr(xml, "BAYERPAT");
    if (!bayer_prop) bayer_prop = strstr(xml, "BayerPattern");
    if (bayer_prop) {
        if (strstr(bayer_prop, "RGGB")) info->bayer_pattern = BAYER_RGGB;
        else if (strstr(bayer_prop, "BGGR")) info->bayer_pattern = BAYER_BGGR;
        else if (strstr(bayer_prop, "GBRG")) info->bayer_pattern = BAYER_GBRG;
        else if (strstr(bayer_prop, "GRBG")) info->bayer_pattern = BAYER_GRBG;
    }

    return 0;
}

// Unshuffle bytes (reverse of XISF byte shuffling)
// Shuffled format: all byte 0s, then all byte 1s, etc.
// Original format: item0[B0,B1,...], item1[B0,B1,...], ...
static void unshuffle_bytes(uint8_t* data, size_t data_size, int item_size) {
    if (item_size <= 1 || data_size == 0) return;

    size_t num_items = data_size / item_size;
    if (num_items == 0) return;

    uint8_t* temp = malloc(data_size);
    if (!temp) return;

    // Copy shuffled data to temp
    memcpy(temp, data, data_size);

    // Unshuffle: each byte position was grouped together
    for (size_t i = 0; i < num_items; i++) {
        for (int b = 0; b < item_size; b++) {
            data[i * item_size + b] = temp[b * num_items + i];
        }
    }

    free(temp);
}

// Decompress data block
static int decompress_block(const uint8_t* compressed, size_t comp_size,
                           uint8_t* output, size_t uncomp_size,
                           XisfCompression codec) {
    switch (codec) {
        case XISF_COMPRESS_NONE:
            if (comp_size != uncomp_size) return -1;
            memcpy(output, compressed, comp_size);
            return 0;

        case XISF_COMPRESS_ZLIB: {
            uLongf dest_len = (uLongf)uncomp_size;
            int ret = uncompress(output, &dest_len, compressed, (uLong)comp_size);
            if (ret != Z_OK) {
                set_xisf_error("zlib decompression failed");
                return -1;
            }
            return 0;
        }

        case XISF_COMPRESS_LZ4:
        case XISF_COMPRESS_LZ4HC: {
            int ret = LZ4_decompress_safe((const char*)compressed, (char*)output,
                                         (int)comp_size, (int)uncomp_size);
            if (ret < 0) {
                set_xisf_error("LZ4 decompression failed");
                return -1;
            }
            return 0;
        }

        case XISF_COMPRESS_ZSTD: {
            size_t ret = ZSTD_decompress(output, uncomp_size, compressed, comp_size);
            if (ZSTD_isError(ret)) {
                char err[256];
                snprintf(err, sizeof(err), "zstd decompression failed: %s", ZSTD_getErrorName(ret));
                set_xisf_error(err);
                return -1;
            }
            return 0;
        }

        default:
            set_xisf_error("Unknown compression codec");
            return -1;
    }
}

// Read raw pixel data from file at given offset
static int read_attachment_data(FILE* f, uint64_t pos, uint64_t size, uint8_t** out_data) {
    *out_data = malloc(size);
    if (!*out_data) {
        set_xisf_error("Out of memory");
        return -1;
    }

    if (fseek(f, (long)pos, SEEK_SET) != 0) {
        set_xisf_error("Failed to seek to attachment position");
        free(*out_data);
        *out_data = NULL;
        return -1;
    }

    size_t read = fread(*out_data, 1, size, f);
    if (read != size) {
        char err[256];
        snprintf(err, sizeof(err), "Failed to read attachment data: expected %llu, got %zu",
                (unsigned long long)size, read);
        set_xisf_error(err);
        free(*out_data);
        *out_data = NULL;
        return -1;
    }

    return 0;
}

// Convert pixel data to float32 for processing (SIMD optimized)
static float* convert_to_float32(const uint8_t* raw_data, const XisfImageInfo* info) {
    size_t num_samples = info->width * info->height * info->channels;
    float* float_data = malloc(num_samples * sizeof(float));
    if (!float_data) {
        set_xisf_error("Out of memory for float conversion");
        return NULL;
    }

    switch (info->sample_format) {
        case XISF_SAMPLE_UINT8: {
            const uint8_t* src = raw_data;
            for (size_t i = 0; i < num_samples; i++) {
                float_data[i] = (float)src[i] * 256.0f;
            }
            break;
        }

        case XISF_SAMPLE_UINT16: {
            const uint16_t* src = (const uint16_t*)raw_data;
            size_t i = 0;
#ifdef HAS_SIMD_SSE2
            // Process 4 uint16 values at a time
            for (; i + 3 < num_samples; i += 4) {
                __m128i v16 = _mm_loadl_epi64((__m128i*)&src[i]);  // Load 4 uint16
                __m128i v32 = _mm_unpacklo_epi16(v16, _mm_setzero_si128());  // Zero-extend to 32-bit
                __m128 vf = _mm_cvtepi32_ps(v32);  // Convert to float
                _mm_storeu_ps(&float_data[i], vf);
            }
#elif defined(HAS_SIMD_NEON)
            for (; i + 3 < num_samples; i += 4) {
                uint16x4_t v16 = vld1_u16(&src[i]);
                uint32x4_t v32 = vmovl_u16(v16);
                float32x4_t vf = vcvtq_f32_u32(v32);
                vst1q_f32(&float_data[i], vf);
            }
#endif
            // Scalar remainder
            for (; i < num_samples; i++) {
                float_data[i] = (float)src[i];
            }
            break;
        }

        case XISF_SAMPLE_UINT32: {
            const uint32_t* src = (const uint32_t*)raw_data;
            for (size_t i = 0; i < num_samples; i++) {
                float_data[i] = (float)(src[i] >> 16);
            }
            break;
        }

        case XISF_SAMPLE_FLOAT32: {
            const float* src = (const float*)raw_data;
            size_t i = 0;
#ifdef HAS_SIMD_SSE2
            __m128 scale = _mm_set1_ps(65535.0f);
            for (; i + 3 < num_samples; i += 4) {
                __m128 v = _mm_loadu_ps(&src[i]);
                v = _mm_mul_ps(v, scale);
                _mm_storeu_ps(&float_data[i], v);
            }
#elif defined(HAS_SIMD_NEON)
            float32x4_t scale = vdupq_n_f32(65535.0f);
            for (; i + 3 < num_samples; i += 4) {
                float32x4_t v = vld1q_f32(&src[i]);
                v = vmulq_f32(v, scale);
                vst1q_f32(&float_data[i], v);
            }
#endif
            for (; i < num_samples; i++) {
                float_data[i] = src[i] * 65535.0f;
            }
            break;
        }

        case XISF_SAMPLE_FLOAT64: {
            const double* src = (const double*)raw_data;
            for (size_t i = 0; i < num_samples; i++) {
                float_data[i] = (float)(src[i] * 65535.0);
            }
            break;
        }
    }

    return float_data;
}

// Convert non-planar RGB to planar format (what our pipeline expects)
static void convert_normal_to_planar(float* data, size_t width, size_t height, size_t channels) {
    if (channels != 3) return;  // Only for RGB

    size_t plane_size = width * height;
    float* temp = malloc(plane_size * 3 * sizeof(float));
    if (!temp) return;

    // Input: RGBRGBRGB... Output: RRR...GGG...BBB...
    for (size_t i = 0; i < plane_size; i++) {
        temp[i] = data[i * 3];                    // R
        temp[i + plane_size] = data[i * 3 + 1];   // G
        temp[i + plane_size * 2] = data[i * 3 + 2]; // B
    }

    memcpy(data, temp, plane_size * 3 * sizeof(float));
    free(temp);
}

int read_xisf_image(const char* xisf_path, FitsMetadata* meta,
                    uint16_t** out_data_u16, float** out_data_f32) {
    FILE* f = fopen(xisf_path, "rb");
    if (!f) {
        set_xisf_error("Failed to open XISF file");
        return -1;
    }

    // Read and validate signature
    uint8_t header[16];
    if (fread(header, 1, 16, f) != 16) {
        set_xisf_error("Failed to read XISF header");
        fclose(f);
        return -1;
    }

    if (memcmp(header, XISF_SIGNATURE, 8) != 0) {
        set_xisf_error("Invalid XISF signature");
        fclose(f);
        return -1;
    }

    // Bytes 8-11: header length (little-endian uint32)
    // Bytes 12-15: reserved
    uint32_t xml_length = header[8] | (header[9] << 8) | (header[10] << 16) | (header[11] << 24);

    if (xml_length == 0 || xml_length > 100 * 1024 * 1024) {  // Sanity check: max 100MB header
        set_xisf_error("Invalid XISF header length");
        fclose(f);
        return -1;
    }

    // Read XML header
    char* xml = malloc(xml_length + 1);
    if (!xml) {
        set_xisf_error("Out of memory for XML header");
        fclose(f);
        return -1;
    }

    if (fread(xml, 1, xml_length, f) != xml_length) {
        set_xisf_error("Failed to read XISF XML header");
        free(xml);
        fclose(f);
        return -1;
    }
    xml[xml_length] = '\0';

    // Parse XML
    XisfImageInfo info;
    if (parse_xisf_xml(xml, xml_length, &info) != 0) {
        free(xml);
        fclose(f);
        return -1;
    }

    // Read pixel data based on location type
    uint8_t* compressed_data = NULL;
    uint8_t* raw_data = NULL;

    if (info.location == XISF_LOCATION_ATTACHMENT) {
        if (read_attachment_data(f, info.attachment_pos, info.block_size, &compressed_data) != 0) {
            free(xml);
            fclose(f);
            return -1;
        }
    } else if (info.location == XISF_LOCATION_EMBEDDED) {
        // Find embedded data in XML (base64 content between > and <)
        const char* data_start = strstr(xml, "<Data");
        if (!data_start) {
            // Try looking for data directly after Image element
            data_start = strstr(xml, ">");
        }
        if (data_start) {
            data_start = strchr(data_start, '>');
            if (data_start) {
                data_start++;
                const char* data_end = strstr(data_start, "</");
                if (data_end) {
                    size_t b64_len = data_end - data_start;
                    size_t decoded_size = base64_decoded_size(b64_len);
                    compressed_data = malloc(decoded_size);
                    if (!compressed_data) {
                        set_xisf_error("Out of memory for embedded data");
                        free(xml);
                        fclose(f);
                        return -1;
                    }
                    size_t actual_len;
                    if (base64_decode(data_start, b64_len, compressed_data, &actual_len) != 0) {
                        set_xisf_error("Base64 decode failed");
                        free(compressed_data);
                        free(xml);
                        fclose(f);
                        return -1;
                    }
                    info.block_size = actual_len;
                }
            }
        }
        if (!compressed_data) {
            set_xisf_error("Failed to find embedded data in XISF");
            free(xml);
            fclose(f);
            return -1;
        }
    } else {
        set_xisf_error("Inline XISF data not yet supported");
        free(xml);
        fclose(f);
        return -1;
    }

    free(xml);
    fclose(f);

    // Decompress if needed
    if (info.compression != XISF_COMPRESS_NONE) {
        raw_data = malloc(info.uncompressed_size);
        if (!raw_data) {
            set_xisf_error("Out of memory for decompression");
            free(compressed_data);
            return -1;
        }

        if (decompress_block(compressed_data, info.block_size,
                            raw_data, info.uncompressed_size, info.compression) != 0) {
            free(compressed_data);
            free(raw_data);
            return -1;
        }
        free(compressed_data);

        // Unshuffle bytes if data was byte-shuffled before compression
        if (info.byte_shuffled && info.shuffle_item_size > 1) {
            unshuffle_bytes(raw_data, info.uncompressed_size, info.shuffle_item_size);
        }
    } else {
        raw_data = compressed_data;
    }

    // Convert to float32 for processing
    float* float_data = convert_to_float32(raw_data, &info);
    free(raw_data);

    if (!float_data) {
        return -1;
    }

    // Convert normal (interleaved) to planar if needed
    if (!info.is_planar && info.channels > 1) {
        convert_normal_to_planar(float_data, info.width, info.height, info.channels);
    }

    // Fill metadata for pipeline compatibility
    meta->width = info.width;
    meta->height = info.height;
    meta->channels = info.channels;
    meta->dtype = DTYPE_FLOAT32;  // We always convert to float32
    meta->bayer_pattern = info.bayer_pattern;
    meta->flip_vertical = info.flip_vertical;

    // Return float data (we've converted everything to float)
    *out_data_u16 = NULL;
    *out_data_f32 = float_data;

    return 0;
}
