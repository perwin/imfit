#include <string.h>   // [PE] for memcpy
#include <stdlib.h>
#include "dream_params.h"

void SetupDreamParams( dream_pars* p, size_t n, const double* init, const string* name,
     	               const int* lock, const double* lo, const double* hi )
{
  p->verboseLevel = 0;
  p->maxEvals = 100000;
  p->numChains = 5;
  p->outputRootname = "";
  p->appendFile = 0;
  p->report_interval = 1;
  p->diagnostics = 0;
  p->burnIn = 0;
  p->recalcLik = 0;
  p->noise = 0.05;            // recommended value in Vrugt+09: 0.05 (alternates: 0.01, 0.1)
  p->bstar_zero = 1e-3;       // recommended value in Vrugt+09: 1.0e-6
  p->collapseOutliers = 1;
  p->gelmanEvals = 5000;
  p->loopSteps = 10;
  p->scaleReductionCrit = 1.01;
  p->deltaMax = 2;
  p->pCR_update = 1;
  p->nCR = 3;                 // recommended value in Vrugt+09: 3
  p->reenterBurnin = 0.2;
  p->fun = NULL;
  p->extraData = NULL;
  p->outputHeaderLines.push_back("# mult_params L burnin gen mult_pCR accept\n");

  p->nvar = n;
  p->nfree = n;
  p->varInit = (double*) calloc(n, sizeof(double));
  p->varLock = (int*) calloc(n, sizeof(int));
  p->varLo = (double*) calloc(n, sizeof(double));
  p->varHi = (double*) calloc(n, sizeof(double));

  memcpy(p->varInit, init,  n*sizeof(double));
  memcpy(p->varLock, lock,  n*sizeof(int));
  memcpy(p->varLo,   lo,    n*sizeof(double));
  memcpy(p->varHi,   hi,    n*sizeof(double));
  for (int i = 0; i < (int)n; ++i) 
    if (lock[i]) 
      --(p->nfree);
  p->deltaMax = (p->nfree - 1)/2;
}



// ---------------------------------------------------------------------------

void SetHeaderDreamParams( dream_pars* p, const vector<string>& headerLine )
{
  p->outputHeaderLines.clear();
  int  nLines = (int)headerLine.size();
  for (int i = 0; i < nLines; i++)
    p->outputHeaderLines.push_back(headerLine[i]);
}



// ---------------------------------------------------------------------------

void FreeVarsDreamParams( dream_pars* p ) {
  free(p->varLo);
  free(p->varHi);
  free(p->varInit);
  free(p->varLock);
  p->nvar = 0;
  p->nfree = 0;
}

