// [intro descriptive comments copied from DESolver.h]
// Differential Evolution Solver Class
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0

// Various modifications by Peter Erwin, 5--16 April 2010
//    We now include a convergence test:
// we sample relative change in best fit statistic (e.g., chi^2) every ten generations; 
// if the relative change is < TOLERANCE (default = 1e-8) for three samples in a row, 
// we declare convergence.

// 19 June 2014: Changed DESolver::Solve to return integer status values instead of bool,
// and also to detect possibly NaN values returned from EnergyFunction [e.g., a signal
// that something went wrong and we should quit] and stop the fitting.

// PROBLEM: DESolver::RandomUniform() uses a constant seed, so it will always generate
// the same sequence of random numbers!


#include <memory.h>
#include <stdio.h>
// Use cmath instead of math.h to avoid GCC-5 problems with C++-11 and isnan()
//#include <math.h>
#include <cmath>
#include <time.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>

using namespace std;

#include "definitions.h"
#include "DESolver.h"
#include "mersenne_twister.h"

// #define DE_CONVERGED   1
// #define DE_MAX_ITERATIONS   5
// #define DE_SIGINT   100
// 

#define Element(a,b,c)  a[b*nDim+c]
#define RowVector(a,b)  (&a[b*nDim])
#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))

const double DEFAULT_TOLERANCE = 1.0e-8;


bool TestConverged( double *relativeDeltas, double ftol );


DESolver::DESolver( int dim, int popSize ) :
          nDim(dim), nPop(popSize),
          generations(0), strategy(stRand1Exp),
          scale(0.7), probability(0.5), trialEnergy(0), bestEnergy(0.0),
          trialSolution(0), bestSolution(0),
          popEnergy(0), population(0), oldValues(0), minBounds(0), maxBounds(0)
{
  trialSolution = new double[nDim];
  bestSolution = new double[nDim];
  popEnergy = new double[nPop];
  population = new double[nPop * nDim];

  // bounds-checking:
  oldValues = new double[nDim];
  minBounds = new double[nDim];
  maxBounds = new double[nDim];
  // tolerance for convergence testing
  tolerance = DEFAULT_TOLERANCE;
}


DESolver::~DESolver( )
{
  if (trialSolution) delete trialSolution;
  if (bestSolution) delete bestSolution;
  if (popEnergy) delete popEnergy;
  if (population) delete population;
  
  if (oldValues) delete oldValues;
  if (minBounds) delete minBounds;
  if (maxBounds) delete maxBounds;

  trialSolution = bestSolution = popEnergy = population = 0;
}


void DESolver::Setup( double *min, double *max, int deStrategy, double diffScale, 
					double crossoverProb, double ftol, unsigned long rngSeed, bool useLHS )
{
  int i;

  strategy = deStrategy;
  scale = diffScale;
  probability = crossoverProb;
  tolerance = ftol;
  
  // PE: seed the (Mersenne Twister) RNG
  if (rngSeed > 0)
    init_genrand(rngSeed);
  else
    init_genrand((unsigned long)time((time_t *)NULL));
  
  CopyVector(minBounds, min);
  CopyVector(maxBounds, max);

  if (useLHS) {
    // Latin hypercube sampling
    printf("   DESolver::Setup -- using Latin hypercube sampling.\n");
    int  sampleOffset;
    double  intervalSize, p;
    // prep Latin hypercube sampling --> sampleIndices
    vector< vector<int> >  sampleIndices;   // [nDim][nPop]
    random_device rd;
    mt19937 g(rd());
    vector<int> singleParamSampleIndices(nPop);
  
    // set up the shuffled indices
    for (int j = 0; j < nDim; j++) {   // iterate over parameters
      for (int i = 0; i < nPop; i++)      // iterate over samples
        singleParamSampleIndices[i] = i;
      shuffle(singleParamSampleIndices.begin(), singleParamSampleIndices.end(), g);
      sampleIndices.push_back(singleParamSampleIndices);
    }

    // generate actual samples
    for (int i = 0; i < nPop; i++) {
      for (int j = 0; j < nDim; j++) {
        sampleOffset = sampleIndices[j][i];
        intervalSize = (max[j] - min[j])/nPop;
        p = genrand_real1();
        Element(population,i,j) = min[j] + (sampleOffset + p)*intervalSize;
      }
      popEnergy[i] = 1.0E20;
    }
  }
  else {
  	// Uniform sampling
    printf("   DESolver::Setup -- using uniform sampling.\n");
 	for (i = 0; i < nPop; i++) {
      for (int j = 0; j < nDim; j++)
        Element(population,i,j) = RandomUniform(min[j], max[j]);
      popEnergy[i] = 1.0E20;
    }
  }

  // old code (uniform Monte Carlo sampling)

  for (i = 0; i < nDim; i++)
    bestSolution[i] = 0.0;
}


// Added by PE
void DESolver::CalcTrialSolution( int candidate )
{
  switch (strategy)
  {
    case stBest1Exp:
      Best1Exp(candidate);
      break;

    case stRand1Exp:
      Rand1Exp(candidate);
      break;

    case stRandToBest1Exp:
      RandToBest1Exp(candidate);
      break;

    case stBest2Exp:
      Best2Exp(candidate);
      break;

    case stRand2Exp:
      Rand2Exp(candidate);
      break;

    case stBest1Bin:
      Best1Bin(candidate);
      break;

    case stRand1Bin:
      Rand1Bin(candidate);
      break;

    case stRandToBest1Bin:
      RandToBest1Bin(candidate);
      break;

    case stBest2Bin:
      Best2Bin(candidate);
      break;

    case stRand2Bin:
      Rand2Bin(candidate);
      break;
  }
}


int DESolver::Solve( int maxGenerations, int verbose )
{
  int generation;
  int candidate;
  bool bAtSolution;
  double  relativeDeltas[3] = {100.0, 100.0, 100.0};
  double  lastBestEnergy;

  bestEnergy = 1.0E20;
  lastBestEnergy = bestEnergy;

  bAtSolution = false;

  for (generation = 0; (generation < maxGenerations) && !bAtSolution; generation++) {
    for (candidate = 0; candidate < nPop; candidate++) {
      // modified by PE
      //(this->*calcTrialSolution)(candidate);
      CalcTrialSolution(candidate);
      // trialSolution now contains a newly generated parameter vector
      // check for out-of-bounds values and generate random values w/in the bounds
      CopyVector(oldValues, RowVector(population, candidate));
      // oldValues is guaranteed to lie between minBounds and maxBounds
      for (int j = 0; j < nDim; j++) {
        if (trialSolution[j] < minBounds[j])
          trialSolution[j] = minBounds[j] + RandomUniform(0.0,1.0)*(oldValues[j] - minBounds[j]);
        if (trialSolution[j] > maxBounds[j])
          trialSolution[j] = maxBounds[j] - RandomUniform(0.0,1.0)*(maxBounds[j] - oldValues[j]);
      }
      
      // Test our newly mutated/bred trial parameter vector
      trialEnergy = EnergyFunction(trialSolution, bAtSolution);

      if (trialEnergy < popEnergy[candidate]) {
        // New low for this candidate
        popEnergy[candidate] = trialEnergy;
        CopyVector(RowVector(population,candidate), trialSolution);

        // Check if all-time low
        if (trialEnergy < bestEnergy) {
          bestEnergy = trialEnergy;
          CopyVector(bestSolution, trialSolution);
        }
      }
    }
    
#ifndef NO_SIGNALS
    if (stopSignal_flag == 1) {
      fprintf(stderr, "\n*** User-requested termination of the fit...\n");
      generations = generation;
      bAtSolution = true;
      return DE_SIGINT;
    }
#endif

    // Debugging printout code added by PE -- print an update every 10 generations
    double  relativeDeltaEnergy = 0.0;
    if ((generation % 10) == 0) {
      if (verbose > 0)
        printf("\nGeneration %4d: bestEnergy = %12.10f", generation, bestEnergy);
      if (generation == 20) {
        relativeDeltaEnergy = fabs(1.0 - lastBestEnergy/bestEnergy);
        relativeDeltas[0] = relativeDeltaEnergy;
        if (verbose > 0)
          printf("   (relative change = %e)", relativeDeltaEnergy);
      }
      else if (generation == 30) {
        relativeDeltaEnergy = fabs(1.0 - lastBestEnergy/bestEnergy);
        relativeDeltas[1] = relativeDeltas[0];
        relativeDeltas[0] = relativeDeltaEnergy;
        if (verbose > 0)
          printf("   (relative change = %e)", relativeDeltaEnergy);
      }
      else if (generation >= 40) {
        relativeDeltaEnergy = fabs(1.0 - lastBestEnergy/bestEnergy);
        relativeDeltas[2] = relativeDeltas[1];
        relativeDeltas[1] = relativeDeltas[0];
        relativeDeltas[0] = relativeDeltaEnergy;
        if (verbose > 0)
          printf("   (relative change = %e)", relativeDeltaEnergy);
        if (TestConverged(relativeDeltas, tolerance)) {
          generations = generation;
          bAtSolution = true;
          return DE_CONVERGED;
        }
      }
      lastBestEnergy = bestEnergy;
    }

	// use "std::isnan" to avoid odd "ambiguity" bug in GCC 4.8.x if you just use "isnan"
    if (std::isnan(bestEnergy)) {
      printf("\n\tcandidate %d, bestEnergy = %f\n", candidate, bestEnergy);
    }

  }
  
  generations = generation;
  return DE_MAX_ITERATIONS;
}


void DESolver::StoreSolution( double *theSolution )
{
  for (int i = 0; i < nDim; i++)
    theSolution[i] = bestSolution[i];
}



void DESolver::Best1Exp( int candidate )
{
  int r1, r2;
  int n;

  SelectSamples(candidate, &r1, &r2);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) {
    trialSolution[n] = bestSolution[n]
              + scale * (Element(population, r1, n)
              - Element(population, r2, n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::Rand1Exp( int candidate )
{
  int r1, r2, r3;
  int n;

  SelectSamples(candidate, &r1, &r2, &r3);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) {
    trialSolution[n] = Element(population, r1, n)
              + scale * (Element(population, r2, n)
              - Element(population, r3, n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::RandToBest1Exp( int candidate )
{
  int r1, r2;
  int n;

  SelectSamples(candidate, &r1, &r2);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) {
    trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
               + scale * (Element(population,r1,n)
               - Element(population,r2,n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::Best2Exp( int candidate )
{
  int r1, r2, r3, r4;
  int n;

  SelectSamples(candidate, &r1, &r2, &r3, &r4);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) {
    trialSolution[n] = bestSolution[n] +
              scale * (Element(population,r1,n)
                    + Element(population,r2,n)
                    - Element(population,r3,n)
                    - Element(population,r4,n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::Rand2Exp( int candidate )
{
  int r1, r2, r3, r4, r5;
  int n;

  SelectSamples(candidate, &r1, &r2, &r3, &r4, &r5);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) {
    trialSolution[n] = Element(population,r1,n)
              + scale * (Element(population,r2,n)
                    + Element(population,r3,n)
                    - Element(population,r4,n)
                    - Element(population,r5,n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::Best1Bin( int candidate )
{
  int r1, r2;
  int n;

  SelectSamples(candidate, &r1, &r2);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; i < nDim; i++) {
    if ((RandomUniform(0.0,1.0) < probability) || (i == (nDim - 1)))
      trialSolution[n] = bestSolution[n]
                + scale * (Element(population,r1,n)
                      - Element(population,r2,n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::Rand1Bin( int candidate )
{
  int r1, r2, r3;
  int n;

  SelectSamples(candidate, &r1, &r2, &r3);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; i < nDim; i++) {
    if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
      trialSolution[n] = Element(population,r1,n)
                + scale * (Element(population,r2,n)
                        - Element(population,r3,n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::RandToBest1Bin( int candidate )
{
  int r1, r2;
  int n;

  SelectSamples(candidate, &r1, &r2);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; i < nDim; i++) {
    if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
      trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
                  + scale * (Element(population,r1,n)
                        - Element(population,r2,n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::Best2Bin( int candidate )
{
  int r1, r2, r3, r4;
  int n;

  SelectSamples(candidate, &r1, &r2, &r3, &r4);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; i < nDim; i++) {
    if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
      trialSolution[n] = bestSolution[n]
                + scale * (Element(population,r1,n)
                      + Element(population,r2,n)
                      - Element(population,r3,n)
                      - Element(population,r4,n));
    n = (n + 1) % nDim;
  }

  return;
}


void DESolver::Rand2Bin( int candidate )
{
  int r1, r2, r3, r4, r5;
  int n;

  SelectSamples(candidate, &r1, &r2, &r3, &r4, &r5);
  n = (int)RandomUniform(0.0, (double)nDim);

  CopyVector(trialSolution, RowVector(population, candidate));
  for (int i = 0; i < nDim; i++) {
    if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
      trialSolution[n] = Element(population,r1,n)
                + scale * (Element(population,r2,n)
                      + Element(population,r3,n)
                      - Element(population,r4,n)
                      - Element(population,r5,n));
    n = (n + 1) % nDim;
  }
}


void DESolver::SelectSamples( int candidate, int *r1, int *r2, int *r3, int *r4, 
								int *r5 )
{
  if (r1) {
    do {
      *r1 = (int)RandomUniform(0.0, (double)nPop);
    } while (*r1 == candidate);
  }

  if (r2) {
    do {
      *r2 = (int)RandomUniform(0.0, (double)nPop);
    } while ((*r2 == candidate) || (*r2 == *r1));
  }

  if (r3) {
    do {
      *r3 = (int)RandomUniform(0.0, (double)nPop);
    } while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
  }

  if (r4) {
    do {
      *r4 = (int)RandomUniform(0.0, (double)nPop);
    } while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
  }

  if (r5) {
    do {
      *r5 = (int)RandomUniform(0.0, (double)nPop);
    } while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)
                          || (*r5 == *r2) || (*r5 == *r1));
  }
}


/// Function added by PE: better random-number-generation function (uses
/// Mersenne Twister and doesn't have a constant seed!)
double DESolver::RandomUniform( double minValue, double maxValue )
{
  double  uniformRand, result;
  
  uniformRand = genrand_real1();   // generates a random number on [0,1]-real-interval
  result = minValue + uniformRand*(maxValue - minValue);
  return result;
}


/// Function added by PE: test for convergence
/// If the last three stored objective-function values (values are stored every 10
/// generations) are all < TOLERANCE, then we decide that we have converged.
bool TestConverged( double *relativeDeltas, double ftol )
{
  if ((relativeDeltas[0] < ftol) && (relativeDeltas[1] < ftol) 
      && (relativeDeltas[2] < ftol))
    return true;
  else
    return false;
}


