#include"./diffraction.h"


int main() {

    Bitmap* slits = importBitmap("slits0.ppm");

    // this does not work for non-square walls.
    Image* pattern = diffract(slits, 1000, 1000, 690, 10);

    exportImage(pattern, "blargh_fraun_brightened.ppm");

    deleteBitmap(slits);
    deleteImage(pattern);
    return 0;
}