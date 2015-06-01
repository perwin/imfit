/** @file
 * \brief Code for printing the results of a fit to stdout or to file.
 *
 *    Utility functions for interpreting and printing results from fits
 *
 */

#ifndef _PRINT_RESULTS_H_
#define _PRINT_RESULTS_H_

#include <string>
#include "mpfit_cpp.h"
#include "model_object.h"
#include "param_struct.h"


/// Code for printing the results of a fit to stdout.
void PrintResults( double *params, mp_result *result, ModelObject *model, 
			int nFreeParameters, mp_par *parameterInfo, int fitStatus );

/// Code for saving the results of a fit to a file.
void SaveParameters( double *params, ModelObject *model, mp_par *parameterInfo, 
          string& outputFilename, string& programName, int argc, char *argv[],
          int nFreeParameters, int whichSolver, int fitStatus );

void SaveParameters2( FILE *file_ptr, double *params, ModelObject *model, mp_par *parameterInfo, 
                    string& programName, int argc, char *argv[], const char *prefix="" );


#endif /* _PRINT_RESULTS_H_ */
