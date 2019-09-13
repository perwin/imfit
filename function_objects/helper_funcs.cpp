#include <math.h>

#include "helper_funcs.h"


const double  A0_M03 = 0.01945;
const double  A1_M03 = -0.8902;
const double  A2_M03 = 10.95;
const double  A3_M03 = -19.67;
const double  A4_M03 = 13.43;


double Calculate_bn( double n )
{
  double  n2 = n*n;
  double  b_n;
  
  if (n > 0.36) {
    // Ciotti & Bertin (1999) approximation, good for all n > 0.36
    b_n = 2*n - 0.333333333333333 + 0.009876543209876543/n
         + 0.0018028610621203215/n2 + 0.00011409410586365319/(n2*n)
         - 7.1510122958919723e-05/(n2*n2);
  } else {
    // MacArthur et al. (2003) approximation for n < 0.36
    b_n = A0_M03 + A1_M03*n + A2_M03*n2 + A3_M03*n2*n + A4_M03*n2*n2;
  }
  return b_n;
}


double CalculateDBEScalingFactor( double h1, double h2, double h3, double r_brk1,
									double r_brk2, double alpha1, double alpha2 )
{
  // n = 2
  double P2a = 1.0 + exp(-alpha1*r_brk1);
  double exp2 = (1.0/alpha1)*(1.0/h1 - 1.0/h2);
  double P2 = pow(P2a, exp2);

  // n = 3
  double P3a = 1.0 + exp(-alpha2*r_brk2);
  double exp3 = (1.0/alpha2)*(1.0/h2 - 1.0/h3);
  double P3 = pow(P3a, exp3);
	
  return 1.0/(P2*P3);
}


double GeneralizedRadius( double deltaX, double deltaY, double cosPA, double sinPA,
							double q, double ellExponent, double invEllExponent )
{
  double  xp, yp_scaled, powerSum;
  // Calculate x,y in component reference frame, and scale y by 1/axis_ratio;
  // convert both to absolute values
  xp = fabs(deltaX*cosPA + deltaY*sinPA);
  yp_scaled = fabs((-deltaX*sinPA + deltaY*cosPA)/q);
  // Compute radius for generalized ellipse
  powerSum = pow(xp, ellExponent) + pow(yp_scaled, ellExponent);
  return pow(powerSum, invEllExponent);
}


double LinearInterp( double r, double r1, double r2, double c01, double c02 )
{
  if (r < r1)
    return c01;
  if (r > r2)
    return c02;
  double m = (c02 - c01) / (r2 - r1);
  return c01 + m*(r - r1);
}
