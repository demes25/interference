#include"./diffraction.h"

// The functions used to calculate intensity are all defined in intensity.h/.c
//
// We may specify which functions we would like to use by defining the relevant 
// macros in diffraction.h
//
// We may also define how much we want the program to 'brighten' the fringes,
// again via macros in the diffraction.h header file.
//
// Units for slit sizes, wall pixel sizes, and wavelengths may be adjusted in intensity.h
#include <time.h>
#include <stdio.h>


int main() {

    Bitmap* slits = importBitmap("images/single_slit/slits.ppm");

    clock_t s = clock();
    Image* pattern = diffract(slits, 1000, 1000, 610, 10);
    clock_t e = clock();

    clock_t dif = e-s;

    println_double((double)dif);

    exportImage(pattern, "images/test.ppm");

    deleteBitmap(slits);
    deleteImage(pattern);
    return 0;
}