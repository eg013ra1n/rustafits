#include "fits_processor.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <string.h>

int save_jpeg(const ProcessedImage* image, const char* output_path, int quality) {
    if (!image || !image->data || !output_path) {
        return -1;
    }

    // Check if output is PNG instead
    size_t len = strlen(output_path);
    if (len > 4 && strcmp(output_path + len - 4, ".png") == 0) {
        stbi_write_png(output_path, image->width, image->height, 3, image->data, image->width * 3);
    } else {
        stbi_write_jpg(output_path, image->width, image->height, 3, image->data, quality);
    }

    return 0;
}
