/** @file
 * \brief Function for estimating memory that imfit will need
 *
 */

#ifndef _ESTIMATE_MEMORY_H_
#define _ESTIMATE_MEMORY_H_

#include <vector>

using namespace std;

#include "psf_oversampling_info.h"


long EstimatePsfOversamplingMemoryUse( vector<PsfOversamplingInfo *> oversamplingInfoVect );

long EstimateMemoryUse( int nData_cols, int nData_rows, int nPSF_cols, int nPSF_rows,
						int nFreeParams, bool levMarFit, bool cashTerms, bool outputResidual,
						bool outputModel );

#endif  // _ESTIMATE_MEMORY_H_
