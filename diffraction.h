/* Demetre Seturidze
 * Arbitrary Slit Interference
 * Diffraction
 * [Header]
 */


#ifndef DIFFRACT
#define DIFFRACT

#include"./intensity.h"

// --- MACROS --- //

// we define our choice of calculation for the intensity.
//
// we have the following choices:
//  
//  DISCRETE -- each pixel on the slit configuration is treated as a point-aperture (no single-slit interference)
//  -   exact_discrete_intensity -- calculates e(ikr) from exact distance
//  -   taylor_discrete_intensity -- approximates distance to second order via Taylor's
//  -   sea_discrete_intensity -- approximates small equal angles, generalizes standard undergrad double-slit formula
//  
//  CONTINUOUS -- various integral approximations are made to account for interference from waves coming from the same slit
//  -   taylor_cerf_intensity -- approximates distance to second order via Taylor's, calculates resulting gaussian integrals using error functions
//  -   fraunhofer_intensity -- uses the fraunhofer diffraction equation (wikipedia), detailed explanation provided with the function definition in intensity.c
#define intensity taylor_cerf_intensity

// we may ask the program to brighten dimmer fringes by applying the 
// transformation u -> pow(u, 1/x). 
// this works because all our final intensity entries u are in [0, 1],
// so dimmer fringes become smoothly brightened for x > 1.
//
// I will define 'brighten' as the x-value noted above.
#define brighten 2.0



// we want to be able to convert the wavelength to a visible color,
// so we can write:
Pixel light_color(double wavelength);

// populates @param wall, @param L meters away from our aperture,
// with the diffraction pattern from the slit configuration in @param slits,
// of monochromatic light with wavelength @param wavelength
Matrix* intensityMatrix(Bitmap* slits, Matrix* wall, double wavelength, double L);


// *** DIFFRACTION *** //

// finds color from wavelength, scales by intensity at each point on wall, 
// returns the resultant image.
Image* diffract(Bitmap* slits, int height, int width, double wavelength, double distance);

#endif