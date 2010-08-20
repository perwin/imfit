/* FILE: print_results.cpp ----------------------------------------- */

#include <stdio.h>
#include <string>

using namespace std;

#include "print_results.h"
#include "mpfit_cpp.h"
#include "param_struct.h"
#include "statistics.h"



/* Local Functions: */
void PrintParam( string& paramName, double paramValue, double paramErr );



// The following is used by PrintResults()
void PrintParam( string& paramName, double paramValue, double paramErr )
{
  if (paramErr == 0.0)
    printf("  %10s = %f\n", paramName.c_str(), paramValue);
  else
    printf("  %10s = %f +/- %f\n", paramName.c_str(), paramValue, paramErr);
}


// This is a function to print the results of a fit.  It's based on code from
// Craig Markwardt's testmpfit.c, but will also accomodate results from a fit
// done with differential evolution (call with result=0 to indicate the latter).
void PrintResults( double *params, double *xact, mp_result *result, ModelObject *model,
									int nFreeParameters, mp_par *parameterInfo, int fitStatus )
{
  int  i;
  int  nValidPixels = model->GetNValidPixels();
  int  nDegreesFreedom = nValidPixels - nFreeParameters;
  string  mpfitMessage;
  double  aic, bic;
  
  if (result == 0) {
    // PrintResult was called with result from Differential Evolution fit, not mpfit
    // Only print results of fit if fitStatus >= 1
    if (fitStatus < 1)
      return;
    double  chiSquared = model->ChiSquared(params);
    printf("  CHI-SQUARE = %lf    (%d DOF)\n", chiSquared, nDegreesFreedom);
    printf("\n");
    aic = AIC_corrected(chiSquared, nFreeParameters, nValidPixels, 1);
    bic = BIC(chiSquared, nFreeParameters, nValidPixels, 1);
    printf("Reduced Chi^2 = %f\n", chiSquared / nDegreesFreedom);
    printf("AIC = %f, BIC = %f\n\n", aic, bic);
    // output the best-fit parameters
    for (i = 0; i < model->GetNParams(); i++) {
      PrintParam(model->GetParameterName(i), params[i] + parameterInfo[i].offset, 0.0);
    }
    return;
  }
  
  // OK, if we got this far, then we're dealing with mpfit output
  InterpretMpfitResult(fitStatus, mpfitMessage);
  printf("*** mpfit status = %d -- %s\n", fitStatus, mpfitMessage.c_str());
    // Only print results of fit if valid fit was achieved
  if ((params == 0) || (result == 0))
    return;
  printf("  CHI-SQUARE = %f    (%d DOF)\n", result->bestnorm, nDegreesFreedom);
  printf("  INITIAL CHI^2 = %f\n", result->orignorm);
  printf("        NPAR = %d\n", result->npar);
  printf("       NFREE = %d\n", result->nfree);
  printf("     NPEGGED = %d\n", result->npegged);
  printf("     NITER = %d\n", result->niter);
  printf("      NFEV = %d\n", result->nfev);
  printf("\n");
  aic = AIC_corrected(result->bestnorm, nFreeParameters, nValidPixels, 1);
  bic = BIC(result->bestnorm, nFreeParameters, nValidPixels, 1);
  printf("Reduced Chi^2 = %f\n", result->bestnorm / nDegreesFreedom);
  printf("AIC = %f, BIC = %f\n\n", aic, bic);
  
  if (xact) {
    for (i = 0; i < result->npar; i++) {
      printf("  P[%d] = %f +/- %f     (ACTUAL %f)\n", 
	     i, params[i], result->xerror[i], xact[i]);
    }
  } else {
    for (i = 0; i < result->npar; i++) {
      PrintParam(model->GetParameterName(i), params[i] + parameterInfo[i].offset, result->xerror[i]);
    }
  }    
}


/* END OF FILE: print_results.cpp ---------------------------------- */
