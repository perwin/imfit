/*! \file
   \brief  Utility functions for string processing, testing file existence, etc.
 */

#ifndef _UTILITIES_MULTIMFIT_H_
#define _UTILITIES_MULTIMFIT_H_

#include <string>
#include <vector>
#include <tuple>

#include "param_struct.h"   // for mp_par structure
#include "model_object.h"
#include "model_object_multimage.h"


void CorrectForImageOffsets( double *paramVect, vector<mp_par> &paramLimits, 
							ModelObjectMultImage *multImageModel );

#endif /* _UTILITIES_MULTIMFIT_H_ */
