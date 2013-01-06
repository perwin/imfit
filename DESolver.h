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

// Minor modifications by Peter Erwin, 5 April 2010

#if !defined(_DESOLVER_H)
#define _DESOLVER_H

#define stBest1Exp			0
#define stRand1Exp			1
#define stRandToBest1Exp	2
#define stBest2Exp			3
#define stRand2Exp			4
#define stBest1Bin			5
#define stRand1Bin			6
#define stRandToBest1Bin	7
#define stBest2Bin			8
#define stRand2Bin			9

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
	~DESolver( void );
	
	// Setup() must be called before solve to set min, max, strategy etc.
	void Setup( double min[], double max[], int deStrategy,
							double diffScale, double crossoverProb, double ftol );

  // CalcTrialSolution is used to determine which strategy to use (added by PE
  // to replace tricky and non-working use of pointers to member functions in
  // original code)
  void CalcTrialSolution( int candidate );
  
	// Solve() returns true if EnergyFunction() returns true.
	// Otherwise it runs maxGenerations generations and returns false.
	virtual bool Solve( int maxGenerations, int verbose=1 );

	// EnergyFunction must be overridden for problem to solve
	// testSolution[] is nDim array for a candidate solution
	// setting bAtSolution = true indicates solution is found
	// and Solve() immediately returns true.
	virtual double EnergyFunction( double testSolution[], bool &bAtSolution ) = 0;
	
	int Dimension(void) { return(nDim); }
	int Population(void) { return(nPop); }

	// Call these functions after Solve() to get results.
	double Energy(void) { return(bestEnergy); }
	
	void StoreSolution( double *theSolution );

	int Generations(void) { return(generations); }

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
