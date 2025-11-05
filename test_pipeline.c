#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "ffi/debayer.h"
#include "ffi/stretch.h"

int main() {
    // Simulate u16 Bayer data (very dark like astronomy image)
    int width = 100;
    int height = 100;
    uint16_t* bayer_data = malloc(width * height * sizeof(uint16_t));

    // Fill with low values (dark image with some bright stars)
    for (int i = 0; i < width * height; i++) {
        if (i % 1000 == 0) {
            bayer_data[i] = 10000;  // Bright star
        } else {
            bayer_data[i] = 100 + (i % 200);  // Dark background
        }
    }

    // Debayer
    int out_width = width / 2;
    int out_height = height / 2;
    float* rgb_data = malloc(out_width * out_height * 3 * sizeof(float));

    super_pixel_rggb_u16(bayer_data, rgb_data, width, height);

    printf("After debayer - first few R channel values:\n");
    for (int i = 0; i < 10; i++) {
        printf("  %f\n", rgb_data[i]);
    }

    // Compute stretch params for R channel
    StretchParams1Channel params;
    compute_stretch_params(rgb_data, out_width, out_height, 65536, &params);

    printf("\nStretch parameters (R channel):\n");
    printf("  Shadows: %f\n", params.shadows);
    printf("  Highlights: %f\n", params.highlights);
    printf("  Midtones: %f\n", params.midtones);
    printf("  Max input: %d\n", params.max_input);

    // Apply stretch
    apply_stretch(rgb_data, out_width, out_height, &params);

    printf("\nAfter stretch - first few R channel values:\n");
    for (int i = 0; i < 10; i++) {
        printf("  %f\n", rgb_data[i]);
    }

    free(bayer_data);
    free(rgb_data);
    return 0;
}
