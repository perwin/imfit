#ifndef __DREAM_PARAMS_H__
#define __DREAM_PARAMS_H__

#include <string>
#include <vector>
using namespace std;

#include "param_struct.h"

typedef double (*LikelihoodFunction)( int chain_id, int gen, const double* state, 
                         const void* pars, bool recalc );

typedef struct t_dream_pars {
  int verboseLevel;                 /* vebose flag   [PE: can be 0, 1, or > 1] */
  int maxEvals;              /* max number of function evaluations */
  int numChains;    
  string outputRootname;             /* output filename */
  int appendFile;            /* continue from previous state */
  int report_interval;       /* report interval for state */
  int diagnostics;           /* report diagnostics at the end of the run */
  int burnIn;                /* number of steps for which to run an adaptive proposal size */
  int recalcLik;             /* recalculate likelihood of previously evaluated states */

  // DREAM variables
  int collapseOutliers;
  int gelmanEvals;            // after burn-in, do Gelman-Rubin convergence test every loopSteps*gelmanEvals interations
  int loopSteps;              // every loopSteps generations, do assorted checks (including convergence)
  int deltaMax;
  int pCR_update;
  int nCR;

  double noise;               // = b in original paper: scaling of gamma = 1 + U(-b,b)
                              // where U(x1,x2) = uniform sampling from specified interval
  double bstar_zero;          // = b^* in original paper: additive offsets = N(0, b^*)
  double scaleReductionCrit;  // declare convergence when Gelman-Rubin scaleReduction
                              // is < scaleReductionCrit for all free variables
  double reenterBurnin;

//               rng->uniform(1, &drand);
//c              e[j] = p->noise*(2.0*drand - 1.0);
//               rng->gaussian(1, &epsilon[j], 0.0, bstar[j]);
//
//  So: e[j] = b * (2*uniform(0,1) - 1)  -- yes, this is e_vector in paper
//      epsilon[j] = N(0, bstar[j])

  int nfree;
  int nvar;
  double* varLo;
  double* varHi;
  double* varInit;
  int* varLock;
  vector<string> parameterNames;

  LikelihoodFunction fun;
  void* extraData;
  
  vector<string> outputHeaderLines;
} dream_pars;


void SetupDreamParams( dream_pars* p, size_t n, const double* init, const string* name,
                    	const int* lock, const double* lo, const double* hi );

// created by PE
void SetHeaderDreamParams( dream_pars* p, const vector<string>& newHeaderLines );

// created by PE
void FreeVarsDreamParams( dream_pars* p );


#endif  // __DREAM_PARAMS_H__
