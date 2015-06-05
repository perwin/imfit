#ifndef _SOLVER_RESULTS_H_
#define _SOLVER_RESULTS_H_

#include "mpfit_cpp.h"

using namespace std;



/// \brief Class for storing useful results related to minimization.
class SolverResults 
{
  public:
    SolverResults( );
    ~SolverResults( );

    void AddMPResults( mp_result& mpResult );

    void SetSolverType( int solverType );
    int GetSolverType( );

    void SetFitStatisticType( int fitStatType );
    int GetFitStatisticType( );
    
    void SetBestfitStatisticValue( double fitStatValue );
    double GetBestfitStatisticValue( );

    void SetNFunctionEvals( int nFunctionEvals );
    int GetNFunctionEvals( );
    
    void SetErrors( double *errors, int nParams );


  private:
    int  whichSolver;
    int  whichFitStatistic;
    double  bestFitValue;
    bool paramSigmasAllocated;
    int nParameters, nFuncEvals;
    double *paramSigmas;
    bool mpResultsPresent;
    mp_result mpResult;
  
};

#endif   // _SOLVER_RESULTS_H_
