/*    Utility functions for interpreting and printing results from fits
 * (useful for 1-D and 2-D fits)
 *
 */

#ifndef _PRINT_RESULTS_H_
#define _PRINT_RESULTS_H_

#include <string>
#include "mpfit_cpp.h"
#include "model_object.h"
#include "param_struct.h"


// Code for printing the results of a fit (either mpfit or differential evolution).
// For a differential-evolution fit, xact and result should be 0; mpStatus will
// be ignored.
void PrintResults( double *params, double *xact, mp_result *result,
					ModelObject *model, int nFreeParameters, mp_par *parameterInfo, int fitStatus );
void SaveParameters( double *params, ModelObject *model, mp_par *parameterInfo, 
          string& outputFilename, string& programName, int argc, char *argv[] );


#endif /* _PRINT_RESULTS_H_ */
