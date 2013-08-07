// Copyright 2010, 2011, 2012, 2013 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "statistics.h"


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

