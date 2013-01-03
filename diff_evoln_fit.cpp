/* FILE: diff_evoln_fit.cpp ---------------------------------------------- */
/*
 * Code for doing Differential Evolution fits; currently only designed for 1-D
 * fitting.
 *
 */

#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DESolver.h"
#include "model_object.h"
#include "param_struct.h"   // for mp_par structure
#include "diff_evoln_fit.h"


#define MAX_DE_GENERATIONS	600



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
  double  lastChiSquared;
  ModelObject  *theModel;
};


double ImfitSolver::EnergyFunction( double *trial, bool &bAtSolution )
{
  double  chiSquared;
  
  chiSquared = theModel->ChiSquared(trial);

  return(chiSquared);
}



// main function called by exterior routines to set up and run the minimization
int DiffEvolnFit( int nParamsTot, double *paramVector, mp_par *parameterLimits, 
                  ModelObject *theModel, double ftol )
{
  ImfitSolver  solver(nParamsTot, 10*nParamsTot, theModel);
  double  *minParamValues;
  double  *maxParamValues;
  int  deStrategy;
  int  maxGenerations;
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
  solver.Setup(minParamValues, maxParamValues, stRandToBest1Exp, F, CR, ftol);

  solver.Solve(maxGenerations);

  solver.StoreSolution(paramVector);

  free(minParamValues);
  free(maxParamValues);
  return 1;
}


/* END OF FILE: diff_evoln_fit.cpp --------------------------------------- */
