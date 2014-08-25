/* FILE: nmsimplex_fit.cpp ----------------------------------------------- */

// Copyright 2014 by Peter Erwin.
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


// Note : the following are the default tolerance values we are currently using
// in mpfitfun.cpp:
//  conf.ftol = 1e-10;   [relative changes in chi^2]
//  conf.xtol = 1e-10;   [relative changes in parameter values]

#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "newlevmar/newlevmar.h"

const int  MAXEVAL_BASE = 10000;
const double  FTOL = 1.0e-8;
const double  XTOL = 1.0e-8;
const int  FUNCS_PER_REPORTING_STEP = 20;
const int  REPORT_STEPS_PER_VERBOSE_OUTPUT = 5;


// Module variables -- used to control user feedback within myfunc
static int  verboseOutput;
static int  funcCount = 0;



// NEW: Object function with interface for levmar code
// p = parameter vector
// e = (unsquared) deviate values = x - hx
// x = data values
// hx = model values
// n = number of data points (pixels)
// m = number of parameters

// FIXME: Needs to actually calculate things properly
//void (*func)( double *p, double *hx, int m, int n, void *adata )
void myfunc_levmar( double *p, double *hx, int m, int n, void *my_func_data )
{
  ModelObject *theModel = (ModelObject *)my_func_data;

  double fitStatistic = theModel->GetFitStatistic(p);

}



// FIXME_EVENTUALLY: currently written to process NLopt return codes
void InterpretResult_newlevmar( int resultValue )
{
  string  description;
  string  returnVal_str;
  ostringstream converter;   // stream used for the conversion
  
  description = "Nelder-Mead Simplex status = ";
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

  printf("%s\n", description.c_str());
}



int NewLevMarFit( int nParamsTot, double *paramVector, mp_par *parameterLimits, 
                  ModelObject *theModel, double ftol, int verbose )
{
  int  maxEvaluations;
  double  finalStatisticVal;
  double  *minParamValues;
  double  *maxParamValues;
  bool  paramLimitsExist = true;
  int  result;
  
  minParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
  maxParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );

  // Check for possible parameter limits
  
  // Specify level of verbosity and start the optimization
  verboseOutput = verbose;
  
  // setup for dlevmar_mle_dif
  int maxIters = 200;
  int nDataVals = theModel->GetNDataValues();
  int fitType = LM_CHISQ_NEYMAN;
  double *dataVals = theModel->GetDataVector();
  result = dlevmar_mle_dif(myfunc_levmar, paramVector, dataVals, nParamsTot, nDataVals,
  				maxIters, NULL, NULL, NULL, NULL, (void *)theModel, fitType);
  
  if (verbose >= 0)
    InterpretResult_newlevmar(result);


  // Dispose of nl_opt object and free arrays:
  free(minParamValues);
  free(maxParamValues);
  return (int)result;
}




/* END OF FILE: nmsimplex_fit.cpp ---------------------------------------- */
