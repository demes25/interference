# interference
A study of various approximation techniques for far-field diffraction patterns with arbitrary slit configurations.

I completed this project for my Optics course in Spring 2025 and recently cleaned and uploaded it to GitHub.

We study here:
- three discrete approximation techniques - where we treat each pixel as a wave point-source located at the upper-left corner (of said pixel)
- two continuous approximation techniques - where we use different methods to approximate the far-field diffraction integral.

This program uses my custom C library clib [https://github.com/demes25/clib] as well as the following external library for computing cerf (complex error function) [https://jugit.fz-juelich.de/mlz/libcerf.git]

## overview

This program acts on a given 'slit map,' which is a binary Portable BitMap (PBM) file. We treat each pixel set to 0 as a 'slit' or aperture for diffraction purposes. When given some specific approximation technique, the program will read the bitmap file and output a Portable PixMap (PPM) file with user-specified dimensions that represents the approximate diffraction pattern from the given slit configuration.

The approximation functions are defined in the intensity.c file.

First, we approximate the total complex wave amplitude at each point on the display wall (functions ending in _waveAt). This is done given the coordinates of the source slit (on the slit map), the coordinates of the target point (on the display wall), the wavelength of incident light, and the separation distance between the apertures and the display wall. (Note: continuous approximations use upper-lower bounds for regions on the slit map, instead of discrete slit coordinates)

The intensity from a slit configuration at a specified point on the display wall is calculated by summing the complex amplitude contribution at the given point from each slit in the configuration and then taking the complex modulus. The input parameters are the slit-map itself, the wavelength of incident light, the coordinates of the desired point on the display wall, and the separation distance between the apertures and the display wall.

There are multiple intensity functions corresponding to various approximations for finding the complex wave amplitude. The desired intensity function may be specified in the diffraction.h header file, by assigning the macro alias 'intensity'. 

We find the diffraction pattern by calculating an intensity matrix. This matrix is populated by scalar values, each entry representing the intensity at the corresponding point on the display wall. We then normalize such that the largest entry in the matrix is 1. Given the wavelength of incident light, we have a conversion function that gives us the corresponding RGB values. We then construct the final diffraction image by scaling the RGB values of our incident light by the corresponding intensity at each point on the display wall (i.e. in the intensity matrix).

We may also 'brighten' the image by taking fractional powers of resultant intensities. Since the intensity matrix is normalized to hold largest value 1, then taking fractional powers (i.e. A^(1/x), for given x) of such entries would smoothly 'brighten' the intensity at each point. This 'brighten' parameter x may be assigned as a macro alias in the diffraction.h header file.

More hints and comments are available in the code files. 

## the problem

In general, the complex wave amplitude with wavenumber k from a single point-slit to a point a distance R away is given by exp(ikR). 
The problem at hand here is approximating the integral Int[A(x, y) exp(ikR) dxdy] over the aperture coordinates x and y, where A(x, y) is some function that indicates the nature of the aperture at the specified (x, y) points. 
The distance R between the slit-point and the display point is a function of the slit coordinates (x, y), the wall coordinates (X, Y), and the separation L between the display wall and the apertures.
By Pythagoras, this distance is given by R = sqrt[(x-X)^2 + (y-Y)^2 + L^2], and it is horrendous to (analytically) integrate something with such a term in the exponential unless we make certain simplifying approximations, which I detail below.

## approximation methods

Here I will briefly detail the approximation methods studied in this project. 

### discrete

Discrete approximations treat each pixel as a point-slit located at the upper-left corner of the pixel. In other words, we treat the aperture function A(x, y), mentioned above, as a collection of delta functions at the upper-left corners of each of the relevant slit pixels. This allows us to simply sum individual discrete values of exp(ikR) evaluated at different R. 

- sea - 'small equal angles', uses a generalized version of the same approximation taught for double-slit interference in standard undergrad physics courses. 
- exact - calculates the exact distance R by using the Pythagorean norm noted.
- taylor - approximates the distance R to second order by assuming that the horizontal separation L is much larger than the vertical distance between the slit points and the display points.

### continuous

Continuous approximations treat each pixel as a full slit and attempt to approximate the integral itself instead of the aperture function. In other words, we take the aperture function A(x, y) to be 1 if (x, y) falls on an open slit, and 0 else. This requires some jumping through hoops.

All of the current continuous approximations assume that the horizontal separation L is much larger than the vertical distance between the slit points and the display points.

- fraunhofer - the Fraunhofer integral approximation is well-substantiated for approximating far-field diffraction. It approximates R to first order, and then approximates the slit coordinates (x, y) to each be much smaller than the display wall coordinates (X, Y). This yields a simple integral of products of exp(ikx) or exp(iky) and acts as a fast continuous approximation.
- taylor + cerf - this approximation is of my own design. It approximates R to second order in the exponential, which leaves a gaussian-like term with a squared exponent, and a higher order term with a quartic exponent. It then approximates the higher-order quartic exponential term to second order via Taylor's. In total this yields a sum of gaussian-weighted polynomials. It then applies integration by parts, along with the complex error function library (cerf) mentioned above, to calculate the relevant integrals of the gaussian-weighted polynomials. 