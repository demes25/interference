/* Demetre Seturidze
 * Arbitrary Slit Interference
 * Intensities
 */

#include"./intensity.h"

// Here I will try to create a program that intakes
// "slits" as a Bitmap (0 if slit, 1 if no slit)
// and outputs an Image that represents the diffraction pattern


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


// ** helpers ** //

// this is a generic intensity calculator for a function of the kind discrete_waveAt(lambda, y, x, Y, X, L)
double _generic_discrete_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L, Complex (*discrete_waveAt)(double, double, double, double, double, double)){
    Complex total_wave = 0;

    for (int i = 0; i < slits -> height; i++)
    for (int j = 0; j < slits -> width; j++)
        // 0's (False) mark the slits
        if (slits -> grid[i][j] == False) {
            double y = (slit_origin_h - i) * slit_unit;
            double x = (slit_origin_w - j) * slit_unit;
            total_wave += discrete_waveAt(lambda, y, x, Y, X, L);
        }
    
    return total_wave * conj(total_wave);
}

// this is a generic intensity calculator for a function of the kind continuous_waveAt(lambda, y_up, y_low, Y, x_up, x_low, X, L)
double _generic_continuous_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L, Complex (*continuous_waveAt)(double, double, double, double, double, double, double, double)){
    Complex total_wave = 0;

    for (int i = 0; i < slits -> height; i++)
    for (int j = 0; j < slits -> width; j++)
        // 0's (False) mark the slits
        if (slits -> grid[i][j] == False) {
            double y_lower = (slit_origin_h - i) * slit_unit;
            double y_upper = y_lower - slit_unit;

            double x_lower = (slit_origin_w - j) * slit_unit;
            double x_upper = x_lower - slit_unit;

            total_wave += continuous_waveAt(lambda, y_lower, y_upper, Y, x_lower, x_upper, X, L);
        }
    
    return total_wave * conj(total_wave);
}

// ---------- //

// the following functions will calculate
// the intensity DISCRETELY.

// i.e. we measure the total intensity from the upper-left 
// corner of each slit pixel.

// we will make several potential functions that do this,
// that make various approximations as to the distance
// between (x, y, 0) and (X, Y, L).


// *** EXACT DISTANCE *** //

// this one does not make ANY approximations as to the distance
// between (x, y, 0) and (X, Y, L)
//
// returns the intensity at point (X, Y) on a wall L meters away,
// given monochromatic light of wavelength lambda and a slit configuration (@param slits)
Complex exact_discrete_waveAt(double lambda, double y, double x, double Y, double X, double L) {
    double dist = sqrt(square(x - X) + square(y - Y) + square(L));
    double s = TAU * dist / lambda; 
    return cexp(I*s);
}

double exact_discrete_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L) {
    return _generic_discrete_intensity(slits, slit_origin_h, slit_origin_w, lambda, Y, X, L, exact_discrete_waveAt);
}


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
Complex taylor_discrete_waveAt(double lambda, double y, double x, double Y, double X, double L) {
    double L2 = 2*L;

    double r2 = (square(x - X) + square(y - Y))/L2;
    double r4 = square(r2)/L2;

    // this is second order in taylor. if we want
    // first order, just take out the r4.
    double s = TAU * (r2 + r4) / lambda; 
    return cexp(I*s);
}

double taylor_discrete_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L) {
    return _generic_discrete_intensity(slits, slit_origin_h, slit_origin_w, lambda, Y, X, L, taylor_discrete_waveAt);
}


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

double sea_discrete_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L) {
    Complex total_wave = 0;

    boolean first = False;
    double first_y, first_x;

    double tanTheta = sqrt(square(Y) + square(X))/L;

    for (int i = 0; i < slits -> height; i++)
    for (int j = 0; j < slits -> width; j++)
        // 0's (False) mark the slits
        if (slits -> grid[i][j] == False) {
            double y = (slit_origin_h - i) * slit_unit;
            double x = (slit_origin_w - j) * slit_unit;

            if (!first) {
                first_y = y; first_x = x;
                total_wave = 1;
                first = True;
            } else {
                double d = sqrt(square(y - first_y) + square(x - first_x));
                double Dphi = TAU * d * tanTheta / lambda;
                total_wave += cexp(I*Dphi);
            }
        }
    
    return total_wave * conj(total_wave);
}




// ---------- //

// here, we will treat the intensity CONTINUOUSLY

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



// this populates the an array of antiderivatives of z^2n exp(-z^2)
void _gauss_antiderivs_upto_even(Complex* result, unsigned int n, Complex z){
    // I[exp(-z^2)] = sqrt(pi)/2  erf(z)
    // i[z exp(-z^2)] = 1/2 (1-exp(-z^2))
    // I[z^n exp(-z^2)] = -z^(n-1)/2 exp(-z^2) + (n-1)/2  I[z^(n-2) exp(-z^2)]   

    // cerf used from the MIT-licensed 'libcerf' library:
    // https://jugit.fz-juelich.de/mlz/libcerf.git
    result[0] = sqrt(PI/2) * cerf(z); 
    
    Complex xp = cexp(-z*z);

    for (int i=1; i <= n; i+= 1){
        double j = 2.0 * i - 1;
        result[i] = 0.5 * (-cpow(z, j) * xp + j*result[i-1]);
    } 
    // we populated this array, where for each [i] we have the 
    // antiderivative of z^2i exp(-z^2)
}

// returns the complex amplitude of the wave on a wall L meters away, at wall coordinates (X, Y);
// by integrating the wave amplitude over the bounds (x_lower, x_upper) x (y_lower, y_upper)
// according to the approximation of the integrand to the second order by taylor,
// and using the complex error function (from libcerf) for gaussian integrals.
Complex taylor_cerf_waveAt(double lambda, double y_lower, double y_upper, double Y, double x_lower, double x_upper, double X, double L){
    // f ~ exp(-z^2)(1 - i/2Lk  z^4 - 1/(8 L^2 k^2)  z^8)) 

    // we write explicitly in terms of xi and chi
    // f ~ exp(-xi^2) exp(-chi^2) (
    //      1 
    //    - i/2Lk 
    //          (xi^4 + 2xi^2 chi^2 + chi^4) 
    //    - 1/(2Lk^2) 
    //          (xi^8 + 4 xi^6 chi^2 + 6 xi^4 chi^4 + 4 xi^2 chi^6 + chi^8))
    //     )

    // we will then need antiderivatives up to order 8, i.e. n = 4.
    // at four points corresponding to two bounds each for two indep variables.
    int n = 4;

    double k = TAU/lambda;
    double _2L = 2*L;

    // we verify that kr^4/(8L^3) << 1
    double condition = k * square(square(x_upper - X) + square(y_upper - Y))/(8*L*L*L);
    if (condition > 0.1){
        println_cstr("Warning: the following value is not << 1");
        println_double(condition);
    }

    Complex z_factor = (1-I) * sqrt(k/(2 * _2L)); // sqrt(-ik/2L) = sqrt(i)* sqrt(k/2L) = [(1-i)/sqrt(2)] sqrt(k/2L)

    // we find integral values by fundl th of calc.
    // - i.e. evaluate the antiderivatives at xi_up, xi_low
    // and then store and work with (xi_up - xi_low)
    Complex xi_up = z_factor * (x_upper - X);
    Complex xi_low = z_factor * (x_lower - X);
    Complex chi_up = z_factor * (y_upper - Y);
    Complex chi_low = z_factor * (y_lower - Y);

    Complex xi_temp[2][n + 1]; 
    Complex chi_temp[2][n + 1];
    // populate temp arrays with antiderivatives
    _gauss_antiderivs_upto_even(xi_temp[1], n, xi_up);
    _gauss_antiderivs_upto_even(xi_temp[0], n, xi_low);
    _gauss_antiderivs_upto_even(chi_temp[1], n, chi_up);
    _gauss_antiderivs_upto_even(chi_temp[0], n, chi_low);

    Complex xi[n+1];
    Complex chi[n+1];
    // subtract antiderivative values at bounds to yield the integral.
    for (int i=0; i <= n ; i+=1){
        xi[i] = xi_temp[1][i] - xi_temp[0][i];
        chi[i] = chi_temp[1][i] - chi_temp[0][i];
    }

    // we write explicitly in terms of xi and chi
    // the complex amplitude f of the light wave is given by approximately
    // f ~ exp(-xi^2) exp(-chi^2) 
    //    - i/2Lk 
    //          ( xi^4 exp(-xi^2) exp(-chi^2) 
    //          + 2 xi^2 exp(-xi^2) exp(-chi^2) chi^2 
    //          + exp(-xi^2) exp(-chi^2) chi^4
    //          ) 
    //    - 1/(8 L^2 k^2) 
    //          ( xi^8 exp(-xi^2) exp(-chi^2) 
    //          + 4 xi^6 exp(-xi^2) exp(-chi^2) chi^2 
    //          + 6 xi^4 exp(-xi^2) exp(-chi^2) chi^4 
    //          + 4 xi^2 exp(-xi^2) exp(-chi^2) chi^6 
    //          + exp(-xi^2) exp(-chi^2) chi^8 
    //          )
    //     )
    //
    // since xi and chi are independent variables, we may
    // integrate the xi/chi dependencies separately, evaluate at 
    // the respective bounds, and then multiply them together.
    //
    // we do this for each additive term, by linearity.
    Complex first_term = xi[0] * chi[0];

    Complex second_term = xi[2] * chi[0] + 2 * xi[1] * chi[1] + xi[0] * chi[2];
    
    Complex third_term = (xi[4] * chi[0] + xi[0] * chi[4]) + 4 * (xi[3] * chi[1] + xi[1] * chi[3]) + 6 * xi[2]*chi[2];

    // we then calculate each order term, multiply with the corresponding overall coefficient,
    // and add together.
    return first_term - second_term * I/(_2L * k) - third_term/(8*L*L*k*k);
    // we return the integral, over the desired source (slit) pixel,
    // of the expression detailed in previous comments.
}

// uses the continuous taylor approximation + cerf for gaussian integrals.
//
// returns the intensity at point (X, Y) on a wall L meters away,
// given monochromatic light of wavelength lambda and a slit configuration (@param slits)
double taylor_cerf_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L) {
    return _generic_continuous_intensity(slits, slit_origin_h, slit_origin_w, lambda, Y, X, L, taylor_cerf_waveAt);
}



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

Complex fraunhofer_waveAt(double lambda, double y_lower, double y_upper, double Y, double x_lower, double x_upper, double X, double L){
    // f ~= Aexp(ikl (X/2 - x)) exp(ikm (Y/2 - y))
    // with l = X/L, m = Y/L
    // 
    // the antiderivative of exp(iks (a - z)) is 
    // i/ks exp(iks(a - z))

    // unless s = 0, in wich case the wave becomes 1
    // and the antiderivative is z.

    double k = TAU/lambda; 

    Complex _x_wave, _y_wave;

    if (X == 0){
        _x_wave = x_upper - x_lower;
    } else { 
        double kl = k*X/L;
        double X_2 = X/2;
        double x_top = X_2 - x_upper;
        double x_bottom = X_2 - x_lower;
        _x_wave = (cexp(I * kl * x_top) - cexp(I*kl*x_bottom))/(kl);
    }

    if (Y == 0) {
        _y_wave = y_upper - y_lower; 
    } else {
        double km = k*Y/L; 
        double Y_2 = Y/2;
        double y_top = Y_2 - y_upper;
        double y_bottom = Y_2 - y_lower;
        _y_wave = (cexp(I * km * y_top) - cexp(I * km * y_bottom))/(km);
    }

    return -(_x_wave * _y_wave);
}

double fraunhofer_intensity(Bitmap* slits, int slit_origin_h, int slit_origin_w, double lambda, double Y, double X, double L){
    return _generic_continuous_intensity(slits, slit_origin_h, slit_origin_w, lambda, Y, X, L, fraunhofer_waveAt);
}

