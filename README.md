# interference
A study of various approximation techniques for far-field diffraction patterns with arbitrary slit configurations.

We study here three discrete approximation techniques - where we treat each pixel as a wave point-source located at the upper-left corner (of said pixel) - and two continuous approximation techniques - where we use different methods to approximate the far-field diffraction integral.
More information about the techniques and formulas is available in the header files intensity.h and diffraction.h.

This program requires the use of clib (c-library) [https://github.com/demes25/clib] as well as the following external library for computing cerf (complex error function) [https://jugit.fz-juelich.de/mlz/libcerf.git]