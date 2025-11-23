/* Demetre Seturidze
 * Arbitrary Slit Interference
 * Diffraction
 */


#include"./diffraction.h"

// we want to be able to convert the wavelength to a visible color,
// so we can write:
Pixel light_color(double wavelength) {
    // conversion algorithm thanks to chatGPT

    byte r = 0, g = 0, b = 0;

    if (wavelength < 400) {
        // UV (380-400 nm): Map to violet/blue range
        r = 0;
        g = (byte)((wavelength - 350) / 50.0 * byte_max);  // Fading from blue to purple
        b = byte_max;
        if (wavelength < 350) {
            b = 128; // More towards blue if closer to 350 nm
        }
    } else if (wavelength >= 400 && wavelength < 450) {
        // Violet (400-450 nm)
        r = (byte)((450 - wavelength) / 50.0 * byte_max);
        b = byte_max;
    } else if (wavelength >= 450 && wavelength < 500) {
        // Blue (450-500 nm)
        g = (byte)((wavelength - 450) / 50.0 * byte_max);
        b = byte_max;
    } else if (wavelength >= 500 && wavelength < 570) {
        // Green (500-570 nm)
        g = byte_max;
        b = (byte)((570 - wavelength) / 70.0 * byte_max);
    } else if (wavelength >= 570 && wavelength < 590) {
        // Yellow (570-590 nm)
        r = (byte)((wavelength - 570) / 20.0 * byte_max);
        g = byte_max;
    } else if (wavelength >= 590 && wavelength < 620) {
        // Orange (590-620 nm)
        r = byte_max;
        g = (byte)((620 - wavelength) / 30.0 * byte_max);
    } else if (wavelength >= 620 && wavelength <= 700) {
        // Red (620-700 nm)
        r = byte_max;
        b = (byte)((wavelength - 620) / 80.0 * byte_max);
    } else if (wavelength > 700 && wavelength <= 750) {
        // IR (700-750 nm): Map to red fading to gray
        r = byte_max;
        g = (byte)((750 - wavelength) / 50.0 * byte_max); // Fade from red to gray
        b = 0;
    } else {
        // Far IR (>750 nm): Fully gray
        r = 128;
        g = 128;
        b = 128;
    }

    return (Pixel){r, g, b};
}

#ifdef brighten
// the following is a helper function for brightening dimmer fringes
// x = 'brighten' is defined as a double greater than 1, so we may brighten by
// applying u -> pow(u, 1/x) to each entry u.
double _brightener(double v){
    return v==0 ? 0 : v==1 ? 1 : pow(v, 1/brighten);    
}

#endif 

// for the intensity matrix, we define our choice of calculation for the intensity:
// populates and returns @param wall
Matrix* intensityMatrix(Bitmap* slits, Matrix* wall, double wavelength, double L) {
    int wall_origin_h = wall -> height/2;
    int wall_origin_w = wall -> width/2;

    int slit_origin_h = slits -> height / 2;
    int slit_origin_w = slits -> width / 2;

    double lambda = wavelength * lambda_unit;

    for (int i = 0; i < wall -> height; i++){
        double Y = (wall_origin_h - i) * wall_unit;
        for (int j = 0; j < wall -> width; j++){
            double X = (wall_origin_w - j) * wall_unit;
            setEntry(wall, i, j, intensity(slits, slit_origin_h, slit_origin_w, lambda, Y, X, L));
        }
    }

    // we will normalize this so that the maximum entry in the matrix is 1.
    double max = 0.0;

    for (nat i = 0; i < wall -> height; i++)
    for (nat j = 0; j < wall -> width; j++)
        if (max < getEntry(wall, i, j)) max = getEntry(wall, i, j);
    
    if (max != 0) bts__scaleMatrix(wall, 1/max, wall);

    // we may include brightness variation if we like.
    // note now that all of our matrix entries are <= 1.
    // so taking fractional powers of each entries should 
    // smoothly brighten dimmer fringes.
    #ifdef brighten

    bts__actOnEntries(wall, _brightener, wall);
    
    #endif 
    
    return wall; 
}


// *** DIFFRACTION *** //

// finds color from wavelength, scales by intensity at each point on wall, 
// returns the resultant image.
Image* diffract(Bitmap* slits, int height, int width, double wavelength, double distance) {
    Matrix* imtx = intensityMatrix(slits, newMatrix(height, width), wavelength, distance);
    
    Image* result = newImage(height, width);
    Pixel color = light_color(wavelength);

    for (nat i = 0; i < height; i++)
    for (nat j = 0; j < width; j++) {
        double local_intensity = getEntry(imtx, i, j);

        result -> grid[i][j].r = (byte)MIN(byte_max, color.r * local_intensity);
        result -> grid[i][j].g = (byte)MIN(byte_max, color.g * local_intensity);
        result -> grid[i][j].b = (byte)MIN(byte_max, color.b * local_intensity);
    }

    deleteMatrix(imtx);
    return result;
}

