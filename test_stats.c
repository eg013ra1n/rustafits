#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "fitsio.h"

int compare_u16(const void* a, const void* b) {
    uint16_t ua = *(const uint16_t*)a;
    uint16_t ub = *(const uint16_t*)b;
    return (ua > ub) - (ua < ub);
}

int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fits_file>\n", argv[0]);
        return 1;
    }

    fitsfile* fptr;
    int status = 0;

    fits_open_file(&fptr, argv[1], READONLY, &status);
    if (status) {
        fprintf(stderr, "Failed to open FITS file\n");
        return 1;
    }

    int naxis;
    long naxes[3];
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 3, naxes, &status);

    size_t npixels = naxes[0] * naxes[1];
    uint16_t* data = malloc(npixels * sizeof(uint16_t));

    int anynul;
    fits_read_img(fptr, TUSHORT, 1, npixels, NULL, data, &anynul, &status);

    printf("Total pixels: %zu\n", npixels);

    // Sort to get percentiles
    uint16_t* sorted = malloc(npixels * sizeof(uint16_t));
    memcpy(sorted, data, npixels * sizeof(uint16_t));
    qsort(sorted, npixels, sizeof(uint16_t), compare_u16);

    printf("Min: %u\n", sorted[0]);
    printf("1%%: %u\n", sorted[npixels / 100]);
    printf("5%%: %u\n", sorted[npixels * 5 / 100]);
    printf("50%% (median): %u\n", sorted[npixels / 2]);
    printf("95%%: %u\n", sorted[npixels * 95 / 100]);
    printf("99%%: %u\n", sorted[npixels * 99 / 100]);
    printf("Max: %u\n", sorted[npixels - 1]);

    free(data);
    free(sorted);
    fits_close_file(fptr, &status);
    return 0;
}
