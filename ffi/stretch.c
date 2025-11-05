/*
    Auto-stretch algorithm based on QuickFits/PixInsight
    Adapted from QuickLook.Plugin.FitsViewer by Siyu Zhang
*/

#include "stretch.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Comparison function for qsort
static int compare_float(const void* a, const void* b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

// Find median using partial sort (nth_element equivalent)
static float median(float* values, int n) {
    int middle = n / 2;
    qsort(values, n, sizeof(float), compare_float);
    return values[middle];
}

void compute_stretch_params(const float* data, int width, int height,
                            int input_range, StretchParams1Channel* params) {
    const int total_pixels = width * height;

    // Sample the data (max 500,000 samples)
    const int max_samples = 500000;
    const int sample_by = total_pixels < max_samples ? 1 : total_pixels / max_samples;
    const int num_samples = total_pixels / sample_by;

    // Allocate sample buffer
    float* samples = (float*)malloc(num_samples * sizeof(float));
    if (!samples) {
        // Fallback to default params
        params->shadows = 0.0f;
        params->highlights = 1.0f;
        params->midtones = 0.5f;
        params->max_input = input_range > 1 ? input_range - 1 : input_range;
        return;
    }

    // Extract samples
    for (int i = 0; i < num_samples; i++) {
        samples[i] = data[i * sample_by];
    }

    // Find median
    float median_sample = median(samples, num_samples);

    // Compute deviations from median
    for (int i = 0; i < num_samples; i++) {
        float diff = samples[i] - median_sample;
        samples[i] = diff >= 0 ? diff : -diff;  // abs
    }

    // Find median deviation
    float med_dev = median(samples, num_samples);
    free(samples);

    // Maximum possible input value
    params->max_input = input_range > 1 ? input_range - 1 : input_range;

    // Normalize to [0, 1] scale
    float normalized_median = median_sample / (float)input_range;
    float MADN = 1.4826f * med_dev / (float)input_range;

    // Determine if we're in upper or lower half
    int upper_half = normalized_median > 0.5f;

    // Compute shadows and highlights
    float shadows, highlights;
    if (upper_half || MADN == 0.0f) {
        shadows = 0.0f;
    } else {
        shadows = fminf(1.0f, fmaxf(0.0f, normalized_median + (-2.8f * MADN)));
    }

    if (!upper_half || MADN == 0.0f) {
        highlights = 1.0f;
    } else {
        highlights = fminf(1.0f, fmaxf(0.0f, normalized_median - (-2.8f * MADN)));
    }

    // Compute midtones
    const float B = 0.25f;
    float X, M;
    if (!upper_half) {
        X = normalized_median - shadows;
        M = B;
    } else {
        X = B;
        M = highlights - normalized_median;
    }

    float midtones;
    if (X == 0.0f) {
        midtones = 0.0f;
    } else if (X == M) {
        midtones = 0.5f;
    } else if (X == 1.0f) {
        midtones = 1.0f;
    } else {
        midtones = ((M - 1.0f) * X) / ((2.0f * M - 1.0f) * X - M);
    }

    // Store parameters
    params->shadows = shadows;
    params->highlights = highlights;
    params->midtones = midtones;
}

void apply_stretch(float* data, int width, int height,
                  const StretchParams1Channel* params) {
    const int total_pixels = width * height;
    const float max_output = 255.0f;
    const int max_input = params->max_input;

    const float midtones = params->midtones;
    const float highlights = params->highlights;
    const float shadows = params->shadows;

    // Precomputed constants
    const float hs_range_factor = (highlights == shadows) ? 1.0f : 1.0f / (highlights - shadows);
    const float native_shadows = shadows * max_input;
    const float native_highlights = highlights * max_input;
    const float k1 = (midtones - 1.0f) * hs_range_factor * max_output / max_input;
    const float k2 = ((2.0f * midtones) - 1.0f) * hs_range_factor / max_input;

    // Apply stretch to each pixel
    for (int i = 0; i < total_pixels; i++) {
        float input = data[i];

        if (input < native_shadows) {
            data[i] = 0.0f;
        } else if (input >= native_highlights) {
            data[i] = max_output;
        } else {
            float input_floored = input - native_shadows;
            data[i] = (input_floored * k1) / (input_floored * k2 - midtones);
        }
    }
}
