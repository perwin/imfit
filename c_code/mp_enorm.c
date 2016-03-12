// Trial versio of mp_enorm(), cleaned up and separated from mpfit.cpp

#include <math.h>

#include "mp_enorm.h"

  /*
   *     **********
   *
   *     function enorm
   *
   *     given an n-vector x, this function calculates the
   *     euclidean norm of x.
   *
   *     the euclidean norm is computed by accumulating the sum of
   *     squares in three different sums. the sums of squares for the
   *     small and large components are scaled so that no overflows
   *     occur. non-destructive underflows are permitted. underflows
   *     and overflows do not occur in the computation of the unscaled
   *     sum of squares for the intermediate components.
   *     the definitions of small, intermediate and large components
   *     depend on two constants, RDWARF and RGIANT. the main
   *     restrictions on these constants are that RDWARF**2 not
   *     underflow and RGIANT**2 not overflow. the constants
   *     given here are suitable for every known computer.
   *
   *     the function statement is
   *
   *  double precision function enorm(n,x)
   *
   *     where
   *
   *  n is a positive integer input variable.
   *
   *  x is an input array of length n.
   *
   *     subprograms called
   *
   *  fortran-supplied ... dabs,dsqrt
   *
   *     argonne national laboratory. minpack project. march 1980.
   *     burton s. garbow, kenneth e. hillstrom, jorge j. more
   *
   *     **********
   *
   * SIMPLE INTERPRETATION: This function returns the square root of the sum
   * of the squares of the elements in vector x.
   */

#define RDWARF   3.834e-20
#define RGIANT   1.304e19


double mp_enorm( int n, double *x )
{
  int  i;
  double  agiant, floatn, xabs;
  double  s1 = 0.0;
  double  s2 = 0.0;
  double  s3 = 0.0;
  double  x1max = 0.0;
  double  x3max = 0.0;
  double  ans, temp;
//  double  rdwarf = RDWARF;
//  double  rgiant = RGIANT;
//  static  double zero = 0.0;
//  static  double one = 1.0;
  
//  s1 = s2 = s3 = 0.0;
//  x1max = x3max = 0.0;
  floatn = n;
  agiant = RGIANT/floatn;
  
  for (i = 0; i < n; i++) {
    xabs = fabs(x[i]);
    
    if ((xabs > RDWARF) && (xabs < agiant)) {
      /* sum for intermediate components. */
      s2 += xabs*xabs;
      continue;
    }
      
    if (xabs > RDWARF) {
      /* sum for large components. */
      if (xabs > x1max) {
        temp = x1max/xabs;
        s1 = 1.0 + s1*temp*temp;
        x1max = xabs;
      }
      else {
        temp = xabs/x1max;
        s1 += temp*temp;
      }
      continue;
    }
    
    /* sum for small components. */
    if (xabs > x3max) {
      temp = x3max/xabs;
      s3 = 1.0 + s3*temp*temp;
      x3max = xabs;
    }
    else {
      if (xabs != 0.0) {
        temp = xabs/x3max;
        s3 += temp*temp;
      }
    }
  }
  
  /* calculation of norm. */
  if (s1 != 0.0) {
    temp = s1 + (s2/x1max)/x1max;
    ans = x1max*sqrt(temp);
    return(ans);
  }
  if (s2 != 0.0) {
    if (s2 >= x3max)
      temp = s2*(1.0 + (x3max/s2)*(x3max*s3));
    else
      temp = x3max*((s2/x3max) + (x3max*s3));
    ans = sqrt(temp);
  }
  else
    ans = x3max*sqrt(s3);
  
  return(ans);
}
