/* Header file for statistics-related functions */

#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#ifdef	__cplusplus
extern "C" {
#endif

double Mean( double *vector, int nVals );

double StandardDeviation( double *vector, int nVals );

void ConfidenceInterval( double *vector, int nVals, double *lower, double *upper );

double AIC_corrected( double logLikelihood, int nParams, long nData, int chiSquareUsed );

double BIC( double logLikelihood, int nParams, long nData, int chiSquareUsed );


#ifdef	__cplusplus
}
#endif

#endif /* _STATISTICS_H_ */
