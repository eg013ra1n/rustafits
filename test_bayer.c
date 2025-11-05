#include <stdio.h>
#include "fitsio.h"

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

    char bayerpat[64] = {0};
    char roworder[64] = {0};

    status = 0;
    fits_read_key(fptr, TSTRING, "BAYERPAT", bayerpat, NULL, &status);
    if (status == 0) {
        printf("BAYERPAT: %s\n", bayerpat);
    } else {
        printf("BAYERPAT: Not found\n");
    }

    status = 0;
    fits_read_key(fptr, TSTRING, "ROWORDER", roworder, NULL, &status);
    if (status == 0) {
        printf("ROWORDER: %s\n", roworder);
    } else {
        printf("ROWORDER: Not found\n");
    }

    fits_close_file(fptr, &status);
    return 0;
}
