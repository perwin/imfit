/* FILE: diff_evoln_fit.cpp ---------------------------------------------- */
/*
 * Code for doing Differential Evolution fits. Implements a subclass of DESolver,
 * specialized for working with ModelObject objects.
 *
 * Main function DiffEvolnFit sets up and runs the fitting process; meant to be 
 * called from other functions (e.g., main() of imfit_main.cpp)
 *
 */

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



#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DESolver.h"
#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "diff_evoln_fit.h"

// "Population" size should be POP_SIZE_PER_PARAMETER * nParametersTot
#define POP_SIZE_PER_PARAMETER  10
#define MAX_DE_GENERATIONS	600

const int  REPORT_STEPS_PER_VERBOSE_OUTPUT = 5;



// Derived DESolver class for our fitting problem

class ImfitSolver : public DESolver
{
public:
  ImfitSolver( int dim, int pop, ModelObject *inputModel ) : DESolver(dim, pop)
  {
    theModel = inputModel;
    count = 0;
  }

  double EnergyFunction( double trial[], bool &bAtSolution );

private:
  int count;
//  double  lastChiSquared;
  ModelObject  *theModel;
};


double ImfitSolver::EnergyFunction( double *trial, bool &bAtSolution )
{
  double  fitStatistic;
  
  fitStatistic = theModel->GetFitStatistic(trial);

  return(fitStatistic);
}



// main function called by exterior routines to set up and run the minimization
int DiffEvolnFit( int nParamsTot, double *paramVector, mp_par *parameterLimits, 
                  ModelObject *theModel, double ftol, int verbose )
{
  ImfitSolver  *solver;
  double  *minParamValues;
  double  *maxParamValues;
  int  deStrategy;
  int  maxGenerations;
  int  nFreeParameters = nParamsTot;
  double  F, CR;   // DE parameters (weight factor (aka "scale"), crossover probability)
  bool  paramLimitsOK = true;
  
  minParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
  maxParamValues = (double *)calloc( (size_t)nParamsTot, sizeof(double) );
  
  // Check for valid parameter limits
  if (parameterLimits == NULL)
    paramLimitsOK = false;
  else {
    for (int i = 0; i < nParamsTot; i++) {
      // user specified a fixed value for this parameter
      if (parameterLimits[i].fixed == 1) {
        minParamValues[i] = paramVector[i];
        maxParamValues[i] = paramVector[i];
        nFreeParameters--;
      }
      else {
        // OK, either we have actual parameter limits, or nothing at all
        if ((parameterLimits[i].limited[0] == 1) && (parameterLimits[i].limited[1] == 1)) {
          // parameter limits for this parameter
          minParamValues[i] = parameterLimits[i].limits[0];
          maxParamValues[i] = parameterLimits[i].limits[1];
        }
        else {
          // oops -- no parameter limits for this parameter!
          paramLimitsOK = false;
        }
      }
    }
  }
  
  if (! paramLimitsOK) {
    printf("\n*** Parameter limits must be supplied for all parameters when using DE!\n");
    printf("Exiting...\n\n");
    free(minParamValues);
    free(maxParamValues);
    return -1;
  }


  // Figure out DE strategy and control parameter values
  deStrategy = stRandToBest1Exp;
  F = 0.85;
  CR = 1.0;
  maxGenerations = MAX_DE_GENERATIONS;
  // Instantiate and set up the DE solver:
//  solver = new ImfitSolver(nParamsTot, POP_SIZE_PER_PARAMETER*nParamsTot, theModel);
  solver = new ImfitSolver(nParamsTot, POP_SIZE_PER_PARAMETER*nFreeParameters, theModel);
  solver->Setup(minParamValues, maxParamValues, stRandToBest1Exp, F, CR, ftol);

  solver->Solve(maxGenerations, verbose);

  solver->StoreSolution(paramVector);

  delete solver;
  free(minParamValues);
  free(maxParamValues);
  return 1;
}


/* END OF FILE: diff_evoln_fit.cpp --------------------------------------- */
