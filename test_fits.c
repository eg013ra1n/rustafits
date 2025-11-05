#include <stdio.h>
#include "fitsio.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fits_file>\n", argv[0]);
        return 1;
    }

    fitsfile* fptr;
    int status = 0;
    int naxis = 0;
    long naxes[3] = {0};
    int bitpix = 0;

    fits_open_file(&fptr, argv[1], READONLY, &status);
    if (status) {
        fprintf(stderr, "Failed to open FITS file (status=%d)\n", status);
        return 1;
    }

    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 3, naxes, &status);
    fits_get_img_type(fptr, &bitpix, &status);

    printf("NAXIS: %d\n", naxis);
    printf("Dimensions: %ld x %ld x %ld\n", naxes[0], naxes[1], naxes[2]);
    printf("BITPIX: %d\n", bitpix);

    switch(bitpix) {
        case BYTE_IMG:      printf("Type: 8-bit unsigned\n"); break;
        case SHORT_IMG:     printf("Type: 16-bit signed\n"); break;
        case LONG_IMG:      printf("Type: 32-bit signed\n"); break;
        case LONGLONG_IMG:  printf("Type: 64-bit signed\n"); break;
        case FLOAT_IMG:     printf("Type: 32-bit float\n"); break;
        case DOUBLE_IMG:    printf("Type: 64-bit double\n"); break;
        case USHORT_IMG:    printf("Type: 16-bit unsigned\n"); break;
        case ULONG_IMG:     printf("Type: 32-bit unsigned\n"); break;
        default:            printf("Type: Unknown\n"); break;
    }

    fits_close_file(fptr, &status);
    return 0;
}
