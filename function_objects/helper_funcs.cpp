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
