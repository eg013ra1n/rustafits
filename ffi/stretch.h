/*
    Auto-stretch algorithm based on QuickFits/PixInsight
    Adapted from QuickLook.Plugin.FitsViewer by Siyu Zhang
*/

#ifndef STRETCH_H
#define STRETCH_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float shadows;
    float highlights;
    float midtones;
    int max_input;
} StretchParams1Channel;

// Compute stretch parameters for one channel
void compute_stretch_params(const float* data, int width, int height,
                            int input_range, StretchParams1Channel* params);

// Apply stretch to one channel (modifies data in-place)
void apply_stretch(float* data, int width, int height,
                  const StretchParams1Channel* params);

#ifdef __cplusplus
}
#endif

#endif /* STRETCH_H */
