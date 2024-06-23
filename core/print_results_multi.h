/** @file
 * \brief Public interface for functions which print the results of a fit 
 *        to stdout or to file -- specialized functions for multimfit.
 *
 *    Utility functions for interpreting and printing results from fits by multimfit.
 *
 */

#ifndef _PRINT_RESULTS_MULTI_H_
#define _PRINT_RESULTS_MULTI_H_

#include <string>
#include "model_object_multimage.h"
#include "definitions_multimage.h"
#include "solver_results.h"


/// Code for saving best-fit model parameters (but not image-description
/// parameters) to a series of files, one per input data image; done via
/// repeated calls to multModel->GetParameterStringsForOneImage
void SaveMultImageParameters( double *params, ModelObjectMultImage *multModel, 
								string& outputFilenameRoot, vector<string>& outputHeader );

/// Code for saving best-fit image-description parameters to a single file
void SaveImageInfoParameters( double *params, ModelObjectMultImage *multModel, 
								vector<ImageInfo>& imageInfoVect, 
								SolverResults& solverResults, string& outputFilename,
								vector<string>& outputHeader );



#endif /* _PRINT_RESULTS_MULTI_H_ */
