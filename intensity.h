/* Demetre Seturidze
 * Arbitrary Slit Interference
 * Intensities
 * [Header]
 */

#ifndef INTENSITY
#define INTENSITY

#include"../library/image.h"
#include"../library/math.h"

// Here I will try to create a program that intakes
// "slits" as a Bitmap (0 if slit, 1 if no slit)
// and outputs an Image that represents the diffraction pattern


#define lambda_unit 1e-9 // nanometer
#define slit_unit 1e-3 // milimeter
#define wall_unit 5e-5 // 50 micrometers

// wavelengths will be in nanometers
static double wavelength;
// lambda will be adjusted to be in meters.
static double lambda;


// the value at a point on the wall (X, Y, L)
// of the wave of light exiting a slit at (x, y, 0) 
// looks like:
// f = Aexp(ik * dist[(X, Y, L), (x, y, 0)])
// where k is 2pi/(wavelength).

// each Bitmap entry will be treated to be mm x mm 
// 
// we will treat each output Pixel to be mm x mm as well.
// so in the exponential, 
// k is in 1/nm; and dist is in mm.
//
// we will take the wall distance L to be in METERS.


// ---------- //

// the following functions will calculate
// the intensity DISCRETELY. 
// MEANING - each pixel is treated as a point-slit at its upper-left corner.
// in the fourier understanding, each pixel is effectively treated as
// a delta function aperture.

// we will make several functions that do this,
// that make various approximations as to the distance
// between (x, y, 0) and (X, Y, L).


// *** EXACT DISTANCE *** //

// this one does not make ANY approximations as to the distance
// between (x, y, 0) and (X, Y, L)

Complex exact_discrete_waveAt(double lambda, double y, double x, double Y, double X, double L);
double exact_discrete_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L);


// *** 2nd ORDER TAYLOR DISTANCE *** //

// this one approximates the distance using a (second-order) taylor approximation.
// namely, it takes the wall distance L to be large such that:
//
// if we let:
// r^2 = (x-X)^2 + (y-Y)^2
// then,
// dist = sqrt(r^2 + L^2) = Lsqrt((r/L)^2 + 1) ~=
// ~= L(1 + (r/L)^2 / 2 + (r/L)^4 / 8) = 
// = L + r^2 / 2L + r^4 / 8L^3
//
// we will remove the L because it represents a constant phase.
// so we will use r^2 / 2L + r^4 / 8L^3
Complex taylor_discrete_waveAt(double lambda, double y, double x, double Y, double X, double L);
double taylor_discrete_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L);



// *** SMALL EQUAL ANGLES *** //

// this one approximates the distance by assuming small angles,
// and that for each point on the wall, the angle it subtends to 
// each point on the slit source is approximately equal.
// namely:
// 
// for X, Y constant, theta is constant.
// and 
// theta ~= sin(theta) ~= tan(theta) = R/L,
// where R^2 = (X^2 + Y^2)
//
// this allows us to explicitly find phase differences.
// we will find these wrt the first slit that we find.
// the phase difference between two slits is given by:
//
// Delta(phi) = k d sin(theta) ~= k d R/L,
// where d is the distance between the slits.
// 
// if the first slit is at x0, y0, then:
//
// d = dist((x0, y0), (x, y))

// this generalizes the standard double-slit formula taught in undergrad
// physics courses.

double sea_discrete_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L);


// ---------- //

// now, we will treat the intensity CONTINUOUSLY

// to do this, we must integrate.
// from each slit, light waves should interfere with each-other!
// for each slit at ij, we integrate from ij to (i+1)(j+1)

// but integrating exactly is difficult, because we have 
// a square root in the exponential.



// *** SECOND ORDER TAYLOR *** //

// We write : f = Aexp(ik sqrt(r^2 + L^2))
//
// We may take a far-field approximation L >> r
// by which
// 
// f = Aexp(ikL sqrt(1 + r^2/L^2)) ~= Aexp(ik(L + r^2/2L + r^4/8L^3)) as per our previous second-order approximation
// 
// we will once again remove the common phase ikL and look at the r-dependent terms
//
// f ~ exp(ik (r^2/2L + r^4/8L^3))
//
// we then allow - xi^2 - chi^2 = ik/2L ((x - X)^2 + (y - Y)^2)
// along with z^2 = xi^2 + chi^2 = -ik/2L r^2
//
// whereby xi = sqrt(-ik/2L) (x - X) and chi = sqrt(-ik/2L) (y - Y)
// 
// so f ~ exp(-z^2 - i/2Lk  z^4)
//
// We will assume that the second term |kr^4/8L^3| << 1, so we may expand the z^4 term 
// as another taylor sum, and exploit the holomorphicity of e^(z^2).
//
// f ~ exp(-z^2)(1 - i/2Lk  z^4 - 1/(2*(2Lk)^2)  z^8))
//
// now by holomorphicity, all these terms have holomorphic antiderivatives 
// which can be found by using the following relations:
// 
// I[exp(-z^2)] = sqrt(pi)/2  erf(z)            -----   [identity]                          
// I[z exp(-z^2)] = 1/2 (1 - exp(-z^2))            -----   [u-substitution]                  
// I[z^n exp(-z^2)] = -z^(n-1)/2 exp(-z^2) + (n-1)/2  I[z^(n-2) exp(-z^2)]            -----   [integ by parts, recursive]
// 
// then we may evaluate integrals by simply plugging things in and subtracting, by fund'l th of calc.


// note here we must calculate the error function of a complex number, which
// is an operation not included in the standard c library. 
// 
// as such, I will be using the libcerf library from MIT
// https://jugit.fz-juelich.de/mlz/libcerf.git
// which I have here pre-compiled into a static .a library.
// the files in ../cerf contain the relevant license text as included
// in the original repository cloned from the above url.
#include"../cerf/cerf.h" 

Complex taylor_continuous_waveAt(double lambda, double x_lower, double x_upper, double X, double y_lower, double y_upper, double Y, double L);
double taylor_continuous_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L);


// *** FRAUNHOFER INTEGRAL *** //

// the Fraunhofer integral gives us the far-field diffraction
// equation by giving an approximation for the complex amplitude of
// the light waves for large L. 

// instead of the above taylor expansion that (only) directly relates r to L,
// the fraunhofer formulation is as follows (wikipedia)

// begin again with f = Aexp(ik sqrt(r^2 + L^2))
// and approximate the exponential as before, but to first order
// and again ignore the common phase ikL

// f ~ Aexp(ikr^2/2L)

// now look at r^2 in terms of the aperture coordinates x,y and the wall coordinates X, Y.
// assuming that |X| >> |x| and |Y| >> |y|, we get the following approximation:

// r^2 = (x-X)^2 + (y-Y)^2 ~= X^2 + Y^2 - 2Xx - 2Yy

// leaving f ~ Aexp(ik[ (X^2 + Y^2)/2L - x X/L - y Y/L ]) = Aexp(ik[ (X/2 - x)X/L + (Y/2 - y)Y/L])
// which we may integrate with great ease, since our apertures are boolean and the exponent is linear in x,y
Complex fraunhofer_waveAt(double lambda, double x_lower, double x_upper, double X, double y_lower, double y_upper, double Y, double L);
double fraunhofer_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L);





#endif