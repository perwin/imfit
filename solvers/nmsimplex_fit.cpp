/* FILE: nmsimplex_fit.cpp ----------------------------------------------- */

// Copyright 2012--2022 by Peter Erwin.
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



// * Return values from nlopt library:
// NLOPT_SUCCESS = 1
// Generic success return value.
// NLOPT_STOPVAL_REACHED = 2
// Optimization stopped because stopval (above) was reached.
// NLOPT_FTOL_REACHED = 3
// Optimization stopped because ftol_rel or ftol_abs (above) was reached.
// NLOPT_XTOL_REACHED = 4
// Optimization stopped because xtol_rel or xtol_abs (above) was reached.
// NLOPT_MAXEVAL_REACHED = 5
// Optimization stopped because maxeval (above) was reached.
// NLOPT_MAXTIME_REACHED = 6
// Optimization stopped because maxtime (above) was reached.
// [edit]
// Error codes (negative return values)
// NLOPT_FAILURE = -1
// Generic failure code.
// NLOPT_INVALID_ARGS = -2
// Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).
// NLOPT_OUT_OF_MEMORY = -3
// Ran out of memory.
// NLOPT_ROUNDOFF_LIMITED = -4
// Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)
// NLOPT_FORCED_STOP = -5
// Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s objective function or constraints.


#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
// Use cmath instead of math.h to avoid GCC-5 problems with C++-11 and isnan()
//#include <math.h>
#include <cmath>

#include <nlopt.h>

#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "nmsimplex_fit.h"
#include "solver_results.h"
#include "utilities_pub.h"
#include "definitions.h"


const int  MAXEVAL_BASE = 10000;
const double  FTOL = 1.0e-8;
const double  XTOL = 1.0e-8;
const int  FUNCS_PER_REPORTING_STEP = 20;
const int  REPORT_STEPS_PER_VERBOSE_OUTPUT = 5;


// Module variables -- used to control user feedback within myfunc_nlopt
static int  verboseOutput;
static int  funcCallCount = 0;
nlopt_opt  optimizer;



/// Objective function: calculates the objective value (ignore gradient calculation)
/// Keep track of how many times this function has been called, and report current
/// chi^2 (or other objective-function value) every 20 calls
/// Note that parameters n and grad are unused, but required by the NLopt interface.
double myfunc_nlopt( unsigned n, const double *x, double *grad, void *my_func_data )
{
  ModelObject *theModel = (ModelObject *)my_func_data;
  // following is a necessary kludge bcs theModel->GetFitStatistic() won't accept const double*
  double  *params = (double *)x;
  double  fitStatistic;
  nlopt_result  junk;
  
  fitStatistic = theModel->GetFitStatistic(params);
  
  // feedback to user
  funcCallCount++;
  if (verboseOutput > 0) {
    if ((funcCallCount % FUNCS_PER_REPORTING_STEP) == 0) {
      printf("\tN-M simplex: function call %d: objective = %f\n", funcCallCount, fitStatistic);
      if ( (verboseOutput > 1) && ((funcCallCount % (REPORT_STEPS_PER_VERBOSE_OUTPUT*FUNCS_PER_REPORTING_STEP)) == 0) ) {
        PrintParametersSimple(theModel, params);
      }
    }
  }
  
  // use "std::isnan" to avoid odd "ambiguity" bug in GCC 4.8.x if you just use "isnan"
  if (std::isnan(fitStatistic)) {
    fprintf(stderr, "\n*** NaN-valued fit statistic detected (N-M optimization)!\n");
    fprintf(stderr, "*** Terminating the fit...\n");
    junk = nlopt_force_stop(optimizer);
  }
#ifndef NO_SIGNALS
  if (stopSignal_flag == 1) {
    fprintf(stderr, "\n*** User-requested termination of the fit...\n");
    junk = nlopt_force_stop(optimizer);
  }
#endif

  return(fitStatistic);
}



void GetInterpretation_NM( const int resultValue, string& outputString )
{
  string  description;
  ostringstream converter;   // stream used for the conversion
  
  description = "Nelder-Mead Simplex: status = ";
  converter << resultValue;      // insert the textual representation of resultValue in the characters in the stream
  description += converter.str();
  
  if (resultValue < 0) {
    description += " -- ERROR:";
    if (resultValue == -1)
      description += " generic (unspecified) failure";
    else if (resultValue == -2)
      description += " invalid arguments!";
    else if (resultValue == -3)
      description += " ran out of memory";
    else if (resultValue == -4)
      description += " roundoff errors limited progress";
    else if (resultValue == -5)
      description += " forced termination called from objective function";
  }
  else if ((resultValue > 0) && (resultValue < 5)) {
    description += " -- SUCCESS:";
    if (resultValue == 1)
      description += " generic (unspecified) success";
    else if (resultValue == 2)
      description += " minimum allowed fit statistic (stopval) reached";
    else if (resultValue == 3)
      description += " ftol_rel or ftol_abs reached";
    else if (resultValue == 4)
      description += " xtol or xtol_abs reached";
  }
  else if (resultValue == 5)
    description += " -- FAILED: reached maximum number of function evaluations";
  else if (resultValue == 6)
    description += " -- FAILED: reached maximum time";

  outputString = description;
}



int NMSimplexFit( const int nParamsTot, double *paramVector, vector<mp_par> parameterLimits, 
                  ModelObject *theModel, const double ftol, const int verbose, 
                  SolverResults *solverResults )
{
  nlopt_result  result;
  int  maxEvaluations;
  double  initialStatisticVal, finalStatisticVal;
  double  *minParamValues;
  double  *maxParamValues;
  
  minParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
  maxParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

  // Extract and store parameter limits, if any
  for (int i = 0; i < nParamsTot; i++) {
    // default state is to have no limits on a parameter
    minParamValues[i] = -HUGE_VAL;
    maxParamValues[i] = HUGE_VAL;
    // check to see if user specified a fixed value for this parameter
    if (parameterLimits[i].fixed == 1) {
      minParamValues[i] = paramVector[i];
      maxParamValues[i] = paramVector[i];
    }
    else if ((parameterLimits[i].limited[0] == 1) && (parameterLimits[i].limited[1] == 1)) {
      // user specified parameter limits for this parameter
      minParamValues[i] = parameterLimits[i].limits[0];
      maxParamValues[i] = parameterLimits[i].limits[1];
    }
  }
  
  // Create an nlopt object, specifying Nelder-Mead Simplex algorithm
  optimizer = nlopt_create(NLOPT_LN_NELDERMEAD, nParamsTot); /* algorithm and dimensionality */
  
  // Specify stopping conditions (desired tolerances, max # function calls)
  // specify relative tolerance (same as ftol in mpfit.cpp)
  nlopt_set_ftol_rel(optimizer, ftol);
  // specify absolute tolerance (same as ftol in mpfit.cpp)
  nlopt_set_ftol_abs(optimizer, ftol);
  // specify relative tolerance for all parameters
  nlopt_set_xtol_rel(optimizer, ftol);
  // maximum number of function calls (MAXEVAL_BASE * total number of parameters)
  maxEvaluations = nParamsTot * MAXEVAL_BASE;
  nlopt_set_maxeval(optimizer, maxEvaluations);
  
  // Set up the optimizer for minimization
  nlopt_set_min_objective(optimizer, myfunc_nlopt, theModel);  
  // Specify parameter boundaries, if they exist
  nlopt_set_lower_bounds(optimizer, minParamValues);
  nlopt_set_upper_bounds(optimizer, maxParamValues);
  
  // record initial fit-statistic value, if we're going to save it
  if (solverResults != NULL)
    initialStatisticVal = theModel->GetFitStatistic(paramVector);
  
  // Specify level of verbosity and start the optimization
  verboseOutput = verbose;
  result = nlopt_optimize(optimizer, paramVector, &finalStatisticVal);
  if (verbose >= 0) {
    string interpretedResult;
    GetInterpretation_NM((int)result, interpretedResult);
    printf("%s\n", interpretedResult.c_str());
  }

  // Store information about the optimization, if SolverResults object was supplied
  if (solverResults != NULL) {
    solverResults->SetSolverType(NMSIMPLEX_SOLVER);
    solverResults->StoreNFunctionEvals(funcCallCount);
    solverResults->StoreBestfitStatisticValue(finalStatisticVal);
    solverResults->StoreInitialStatisticValue(initialStatisticVal);
  }

  // Dispose of nl_opt object and free arrays:
  nlopt_destroy(optimizer);
  free(minParamValues);
  free(maxParamValues);
  return (int)result;
}



/* END OF FILE: nmsimplex_fit.cpp ---------------------------------------- */
