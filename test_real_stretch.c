#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "ffi/stretch.h"

int main() {
    // Load some sample RGB channel data (simulating dark astronomy image)
    // Let's create a test with similar values to real FITS
    int width = 100;
    int height = 100;
    float* data = malloc(width * height * sizeof(float));

    // Simulate very dark image with a few bright spots
    for (int i = 0; i < width * height; i++) {
        if (i < 10) {
            data[i] = 5000.0f + i * 100;  // Bright pixels
        } else if (i < 100) {
            data[i] = 100.0f + i;  // Medium pixels
        } else {
            data[i] = 50.0f + (i % 100);  // Dark background
        }
    }

    printf("Sample input values:\n");
    printf("  Bright: %.1f, %.1f, %.1f\n", data[0], data[5], data[9]);
    printf("  Medium: %.1f, %.1f, %.1f\n", data[10], data[50], data[99]);
    printf("  Dark: %.1f, %.1f, %.1f\n", data[100], data[500], data[5000]);

    // Compute stretch params
    StretchParams1Channel params;
    compute_stretch_params(data, width, height, 65536, &params);

    printf("\nComputed stretch parameters:\n");
    printf("  Shadows: %f (native: %.1f)\n", params.shadows, params.shadows * params.max_input);
    printf("  Highlights: %f (native: %.1f)\n", params.highlights, params.highlights * params.max_input);
    printf("  Midtones: %f\n", params.midtones);
    printf("  Max input: %d\n", params.max_input);

    // Apply stretch
    apply_stretch(data, width, height, &params);

    printf("\nSample output values (after stretch, 0-255 range):\n");
    printf("  Bright: %.1f, %.1f, %.1f\n", data[0], data[5], data[9]);
    printf("  Medium: %.1f, %.1f, %.1f\n", data[10], data[50], data[99]);
    printf("  Dark: %.1f, %.1f, %.1f\n", data[100], data[500], data[5000]);

    free(data);
    return 0;
}
