#include "ffi/stretch.h"
#include <stdio.h>

int main() {
    // Test with some sample data
    float test_data[100];
    for (int i = 0; i < 100; i++) {
        test_data[i] = i * 100.0f;  // 0 to 9900
    }

    StretchParams1Channel params;
    compute_stretch_params(test_data, 10, 10, 65536, &params);

    printf("Test stretch parameters:\n");
    printf("  Shadows: %f\n", params.shadows);
    printf("  Highlights: %f\n", params.highlights);
    printf("  Midtones: %f\n", params.midtones);
    printf("  Max input: %d\n", params.max_input);

    return 0;
}
