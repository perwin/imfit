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

// Minor modifications by Peter Erwin, 5 April 2010; 19 June 2014

#ifndef _DESOLVER_H
#define _DESOLVER_H

// #if !defined(_DESOLVER_H)
// #define _DESOLVER_H

const int stBest1Exp       =    0;
const int stRand1Exp       =    1;
const int stRandToBest1Exp =    2;
const int stBest2Exp       =    3;
const int stRand2Exp       =    4;
const int stBest1Bin       =    5;
const int stRand1Bin       =    6;
const int stRandToBest1Bin =    7;
const int stBest2Bin       =    8;
const int stRand2Bin       =    9;

class DESolver;

// this defines a type called "StrategyFunction" which is a pointer to
// a member function of DESolver, which takes an int and returns void
// Currently commented out bcs compilation errors resulted when trying to
// use it.
//typedef void (DESolver::*StrategyFunction)(int);

class DESolver
{
public:
  DESolver( int dim, int popSize );
  ~DESolver( );
	
  /// Setup() must be called before Solve to set min, max, strategy etc.
  void Setup( double min[], double max[], int deStrategy,
							double diffScale, double crossoverProb, double ftol );

  /// CalcTrialSolution is used to determine which strategy to use (added by PE
  /// to replace tricky and non-working use of pointers to member functions in
  /// original code)
  void CalcTrialSolution( int candidate );
  
  virtual int Solve( int maxGenerations, int verbose=1 );

  // EnergyFunction must be overridden for problem to solve
  // testSolution[] is nDim array for a candidate solution
  // setting bAtSolution = true indicates solution is found
  // and Solve() immediately returns true.
  virtual double EnergyFunction( double testSolution[], bool &bAtSolution ) = 0;
	
  int Dimension( ) { return(nDim); }

  int Population( ) { return(nPop); }

  /// Call after Solve() to get results.
  double Energy( ) { return(bestEnergy); }
	
  /// Call after Solve() to get results.
  void StoreSolution( double *theSolution );

  /// Call after Solve() to get results.
  int Generations( ) { return(generations); }

protected:
  void SelectSamples( int candidate, int *r1, int *r2=0, int *r3=0, 
												int *r4=0, int *r5=0 );
  double RandomUniform( double min, double max );

  int nDim;
  int nPop;
  int generations;

  int strategy;
//	StrategyFunction calcTrialSolution;
  double scale;
  double probability;

  double trialEnergy;
  double bestEnergy;

  double *trialSolution;
  double *bestSolution;
  double *popEnergy;
  double *population;

  // added by PE for bounds-checking
  double *oldValues;
  double *minBounds;
  double *maxBounds;
  // added by PE for user specification of fractional tolerance (for convergence test)
  double  tolerance;
	
	// added by PE for debugging purposes
  double  lastBestEnergy;

private:
  void Best1Exp(int candidate);
  void Rand1Exp(int candidate);
  void RandToBest1Exp(int candidate);
  void Best2Exp(int candidate);
  void Rand2Exp(int candidate);
  void Best1Bin(int candidate);
  void Rand1Bin(int candidate);
  void RandToBest1Bin(int candidate);
  void Best2Bin(int candidate);
  void Rand2Bin(int candidate);
};

#endif // _DESOLVER_H
