#ifndef __DREAM_H__
#define __DREAM_H__

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>

#include <rng/RngStream.h>
#include "array.h"
#include "dream_params.h"


// Return values
const int  DREAM_EXIT_NO_OUTPUT_FILES = -10;
const int  DREAM_EXIT_ERROR = -1;
const int  DREAM_EXIT_CONVERGENCE = 0;
const int  DREAM_EXIT_MAX_ITERATIONS = 1;


int dream_restore_state( const dream_pars* p, Array3D<double>& state, 
						Array2D<double>& lik, vector<double>& pCR, int& inBurnIn );

void dream_initialize( const dream_pars* p, rng::RngStream* rng, 
						Array2DView<double>& state, ArrayView<double>& lik );

int dream( const dream_pars* p, rng::RngStream* rng );

void check_outliers( int t, Array2D<double>& lik, vector<double>& meanlik,
                    vector<bool>& outliers );

void gen_CR( rng::RngStream* rng, const vector<double>& pCR, 
            Array2D<int>& CRm, vector<unsigned>& L );

void gelman_rubin( Array3D<double>& state, vector<double>& scaleReduction, 
                  const int* lockVar, int first = 0, int numIter = -1, 
                  int adjustDF = 0 );

#endif   // __DREAM_H__
