#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "statistics.h"
#include "mersenne_twister.h"


/* lower and upper bounds of 68.3% confidence interval: */
#define ONESIGMA_LOWER   0.1585
#define ONESIGMA_UPPER   0.8415



int SortComp(const void *x, const void *y );



/* ---------------- Mean ----------------------------------------------- */

double Mean( double *vector, int nVals )
{
  int  n;
  double  sum = 0.0;
  
  for (n = 0; n < nVals; n++)
    sum += vector[n];
    
  return sum/nVals;
}




/* ---------------- StandardDeviation ---------------------------------- */

double StandardDeviation( double *vector, int nVals )
{
  int  n;
  double  mean = 0.0;
  double  sum = 0.0;
  double  diff;
  
  for (n = 0; n < nVals; n++)
    mean += vector[n];
  mean = mean / nVals;
  for (n = 0; n < nVals; n++) {
    diff = vector[n] - mean;
    sum += diff*diff;
  }
  if (sum >= 0.0)
    return sqrt( sum / (nVals - 1) );
  else
    return 0.0;
}




/* ---------------- ConfidenceInterval --------------------------------- */
/*    Returns lower and upper bounds of a confidence interval for a range
 * of values.  For now, we assume a 68.3% confidence interval (which for
 * a Gaussian distribution is the +/- 1-sigma interval), though in the future
 * this could be an input parameter.
 *    Note that the input vector is sorted in place!
 */
void ConfidenceInterval( double *vector, int nVals, double *lower, double *upper )
{
  int  lower_ind, upper_ind;
  
  lower_ind = round(ONESIGMA_LOWER * nVals) - 1;
  upper_ind = round(ONESIGMA_UPPER * nVals);
  
  /* Sort the vector into increasing order: */
  qsort((void*)vector, (size_t)nVals, sizeof(vector[0]), SortComp);
  *lower = vector[lower_ind];
  *upper = vector[upper_ind];
}




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



/* ---------------- SortComp ------------------------------------------- */
/* Callback function used by qsort() in ConfidenceInterval() to sort a
 * vector of doubles
 */
int SortComp(const void *x, const void *y )
{
  const double  *xx = (const double *)x;
  const double  *yy = (const double *)y;
  
  if (*xx > *yy)
    return 1;
  else {
    if (*xx == *yy)
      return 0;
    else
      return -1;
  }
}



/* ---------------- AIC_corrected -------------------------------------- */
/* Calculate the bias-corrected Akaike Information Criterion for a model fit
 * to data, given the ln(likelihood) of the best-fit model, the number of model
 * parameters nParams, and the number of data points nData (the latter is used
 * to correct the 2*nParams part of AIC for small sample size).
 *
 * If chiSquareUsed is nonzero, then the input logLikelihood is assumed to be
 * the chi^2 value of the fit, and that is used for -2 ln(likelihood)
 * in the calculation.
 *
 * Formula from Burnham & Anderson, Model selection and multimodel inference: 
 * a practical information-theoretic approach (2002), p.66.
 */
double AIC_corrected( double logLikelihood, int nParams, long nData, int chiSquareUsed )
{
  double  twok, aic, correctionTerm;
  
  twok = 2.0*nParams;
  if ( chiSquareUsed )  // user passed chi^2 value as "logLikelihood"
    aic = logLikelihood + twok;
  else  // "logLikelihood" really is ln(likelihood)
    aic = -2.0*logLikelihood + twok;

  correctionTerm = 2.0*nParams*(nParams + 1.0) / (nData - nParams - 1.0);
  return (aic + correctionTerm);
}



/* ---------------- BIC ------------------------------------------------ */
/* Calculate the Bayesian Information Criterion for a model fit to data,
 * given the ln(likelihood) of the best-fit model, the number of model 
 * parameters nParams, and the number of data points nData.
 *
 * If chiSquareUsed is nonzero, then the input logLikelihood is assumed to be
 * the chi^2 value of the fit, and that is used for -2 ln(likelihood)
 * in the calculation.
 */
double BIC( double logLikelihood, int nParams, long nData, int chiSquareUsed )
{
  double  minustwo_logLikelihood, bic;
  
  if ( chiSquareUsed )  // user passed chi^2 value as "logLikelihood"
    minustwo_logLikelihood = logLikelihood;
  else  // "logLikelihood" really is ln(likelihood)
    minustwo_logLikelihood = -2.0*logLikelihood;

  //printf("nParams = %d, nData = %ld\n", nParams, nData);
  bic = minustwo_logLikelihood + nParams*log(nData);
  return bic;
}

