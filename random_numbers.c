#include <math.h>
#include "mersenne_twister.h"
#include "random_numbers.h"


/* ---------------- GaussianRand --------------------------------------- */
/* Based on gasdev() from Numerical Recipes in C, 2nd ed. (p.289).
 * Note that we are currently using genrand_real3(), the (0,1) rng from the
 * Mersenne Twister family (mersenne_twister.c). */

double GaussianRand( void )
{
  static int  iset = 0;
  static double  gset;
  double  fac, rsq, v1, v2;
  
  if (iset == 0) {
    do {
      v1 = 2.0*genrand_real3() - 1.0;
      v2 = 2.0*genrand_real3() - 1.0;
      rsq = v1*v1 + v2*v2;
    } while ( (rsq >= 1.0) || (rsq == 0.0) );
    fac = sqrt(-2.0*log(rsq)/rsq);
    
    /* Box-Muller transformation to get two normal deviates; return one and
     * save the other for the next call. */
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}
