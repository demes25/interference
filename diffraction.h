/* Demetre Seturidze
 * Arbitrary Slit Interference
 * Diffraction
 * [Header]
 */


#ifndef DIFFRACT
#define DIFFRACT

#include"./intensity.h"

// we want to be able to convert the wavelength to a visible color,
// so we can write:
Pixel light_color(double wavelength);

// for the intensity matrix, we define our choice of calculation for the intensity:
#define intensity fraunhofer_intensity

// we may define this to brighten dimmer fringes by applying fractional powers
// (say sqrt). this works because all our final intensity entries are in [0, 1],
// so dimmer fringes become smoothly brightened.
//
// we will apply mtx[ij] -> pow(mtx[ij], 1/x)
// where x is the value of the macro 'brighten'
#define brighten 2.0

Matrix* intensityMatrix(Bitmap* slits, Matrix* wall, double wavelength, double L);


// *** DIFFRACTION *** //

// finds color from wavelength, scales by intensity at each point on wall, 
// returns the resultant image.
Image* diffract(Bitmap* slits, int height, int width, double wavelength, double distance);

#endif