/* FILE: func1d_n1543majmin2.cpp --------------------------------------- */
/* VERSION 0.1
 *
 *   Function object class for fitting major + minor-axis profile data for NGC 1543.
 * Includes extra central Gaussian component (for, e.g., nuclear star cluster).
 *
 * Flux-related parameters are in surface-brightness (mag/arcsec^2), but output
 * is flux (calling function -- e.g., ModelObject1D::CreateModelImage -- will
 * convert back to magnitudes for comparison with data
 ;    Parameters: p[0] = n, p[1] = mu_e, p[2] = r_e_maj, p[3] = r_e_min, p[4] = mu_0_bar,
;                p[5] = h_bar1, p[6] = h_bar2, p[7] = r_b, p[8] = alpha,
;                p[9] = mu_0_disk, p[10] = h_disk_maj, p[11] = h_disk_min,
;                p[12] = sigma_bar, and p[13] = x_major_max

 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x.
 *      GetValue() then completes the calculation, using the actual value
 *      of x, and returns the result.
 *
 *   MODIFICATION HISTORY:
 */


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func1d_n1543majmin2.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 15;
const char  PARAM_LABELS[][20] = {"n", "mu_e", "r_e_maj", "r_e_min", "mu0_bar", "h_bar1", 
	"h_bar2", "r_b", "alpha", "mu0_disk", "h_disk_maj", "h_disk_min", "sigma_bar",
	"mu0_nsc", "sigma_nsc"};
const char FUNCTION_NAME[] = "n1543majmin2-1D function";
#define CLASS_SHORT_NAME  "n1543majmin2-1D"

const char n1543MajMin21D::className[] = CLASS_SHORT_NAME;


/* ---------------- CONSTRUCTOR ---------------------------------------- */

n1543MajMin21D::n1543MajMin21D( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = CLASS_SHORT_NAME;

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void n1543MajMin21D::Setup( double params[], int offsetIndex, double xc )
{
  x0 = xc;
  n = params[0 + offsetIndex];
  mu_e = params[1 + offsetIndex];
  r_e_maj = params[2 + offsetIndex];
  r_e_min = params[3 + offsetIndex];
  mu_0_bar = params[4 + offsetIndex ];
  h_bar1 = params[5 + offsetIndex ];
  h_bar2 = params[6 + offsetIndex ];
  r_b = params[7 + offsetIndex ];
  alpha = params[8 + offsetIndex ];
  mu_0_disk = params[9 + offsetIndex ];
  h_disk_maj = params[10 + offsetIndex ];
  h_disk_min = params[11 + offsetIndex ];
  sigma_bar = params[12 + offsetIndex ];
  mu_0_nsc = params[13 + offsetIndex ];
  sigma_nsc = params[14 + offsetIndex ];
  
  // pre-compute useful things for this round of invoking the function
  
  // Sersic stuff
  I_e = pow(10.0, 0.4*(ZP - mu_e));
  n2 = n*n;
  /* The following approximation for b_n is good for all
   * n > 0.36 */
  bn = 2*n - 0.333333333333333 + 0.009876543209876543/n
       + 0.0018028610621203215/n2 + 0.00011409410586365319/(n2*n)
       - 7.1510122958919723e-05/(n2*n2);
  invn = 1.0 / n;
  
  // Exp + Broken-exp. stuff
  I_0_disk = pow(10.0, 0.4*(ZP - mu_0_disk));
  I_0_bar = pow(10.0, 0.4*(ZP - mu_0_bar));
  exponent = (1.0/alpha) * (1.0/h_bar1 - 1.0/h_bar2);
  // Calculate S [note that this can cause floating *underflow*, but that's OK]:
  S = pow( (1.0 + exp(-alpha*r_b)), (-exponent) );
  // Gaussian for NSC
  I_0_nsc = pow(10.0, 0.4*(ZP - mu_0_nsc));
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double n1543MajMin21D::GetValue( double x )
{
  double  r = fabs(x - x0);
  double  I = 0.0;
  
  // We handle the two different pieces of the profile by assigning *negative*
  // values to the *minor*-axis radii
  
  if (x < 0.0) {
    // Minor-axis fitting
    I += Sersic(r, n, I_e, r_e_min);
    I += Exponential(r, I_0_disk, h_disk_min);
    I += Gaussian(r, I_0_bar, sigma_bar);
  }
  else {
    // Major-axis fitting
    I += Sersic(r, n, I_e, r_e_maj);
    I += Exponential(r, I_0_disk, h_disk_maj);
    // broken-exponential code
    if ( alpha*(r - r_b) > 100.0 ) {
      // Outer-exponential approximation:
      I += I_0_bar * S * exp(r_b/h_bar2 - r_b/h_bar1 - r/h_bar2);
    } else {
      // no danger of overflow in exponentiation, so use fully correct calculation:
      I += I_0_bar * S * exp(-r/h_bar1) * pow( 1.0 + exp(alpha*(r - r_b)), exponent );
  }
  // Finally, add nuclear star cluster Gaussin
  I += Gaussian(r, I_0_nsc, sigma_nsc);
//  printf("In GetValue: x = %g, I_0_bar = %g, h = %g\n", x, I_0_bar, h);
  }
  return (I);
}




/* Extra auxiliary functions: */

/* ---------------- PRIVATE METHOD: Gaussian ---------------------------- */

double n1543MajMin21D::Gaussian( double x, double I_0, double sigma )
{
  // We assume that x >= 0
  double  scaledDeltaR = x / sigma;
  double  expon = -(scaledDeltaR * scaledDeltaR)/2.0;
  return I_0 * exp( expon );
}


/* ---------------- PRIVATE METHOD: Exponential ------------------------- */

double n1543MajMin21D::Exponential( double x, double I_0, double h )
{
  // We assume that x >= 0
  return (I_0 * exp(-x/h));
}


/* ---------------- PRIVATE METHOD: Sersic ----------------------------- */

double n1543MajMin21D::Sersic( double x, double nn, double Ie, double re )
{
  // We assume that x >= 0
  double  I = Ie * exp( -bn * (pow((x/re), invn) - 1.0) );
  return (I);
}



/* END OF FILE: func1d_n1543majmin2.cpp -------------------------------- */
