// Miscellaneous functions used by some of the FunctionObject subclasses,
// extracted and placed here to avoid code duplication.

#ifndef _HELPER_H_
#define _HELPER_H_


// Sersic-related

/// Calculate the b_n parameter for a Sersic function
double Calculate_bn( double n );


/// Calculate scaling factor for double-broken-exponential ("DBE")
double CalculateDBEScalingFactor( double h1, double h2, double h3, double r_brk1,
									double r_brk2, double alpha1, double alpha2 );


// Generalized ellipse shapes

/// Calculate equivalent radius for a generalized ellipse
///   deltaX = x - x0, deltaY = y - y0
///   cosPA, sinPA = cosine and sine of PA_rad, where
///      PA_rad = (PA + 90.0) converted to radians
///   q = axis ratio b/a of ellipse
///   ellExponent = c0 + 2, invEllExponent = 1/ellExponent
double GeneralizedRadius( double deltaX, double deltaY, double cosPA, double sinPA,
							double q, double ellExponent, double invEllExponent );



// Experimental functions for interpolating c0 values

double LinearInterp( double r, double r1, double r2, double c01, double c02 );


// Hyperbolic-tangent truncation function
double HyperbolicTangentTrunc( double r, double r_brk, double B );


#endif  // _HELPER_H_
