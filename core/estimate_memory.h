/** @file
 * \brief Function for estimating memory that imfit will need
 *
 */

#ifndef _ESTIMATE_MEMORY_H_
#define _ESTIMATE_MEMORY_H_

long EstimateMemoryUse( int nData_cols, int nData_rows, int nPSF_cols, int nPSF_rows,
						int nFreeParams, bool levMarFit, bool cashTerms, bool outputResidual,
						bool outputModel, int nPSF_osamp_cols=0, int nPSF_osamp_rows=0,
						int deltaX_osamp=0, int deltaY_osamp=0, int oversampleScale=0 );

#endif  // _ESTIMATE_MEMORY_H_
