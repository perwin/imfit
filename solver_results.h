#ifndef _SOLVER_RESULTS_H_
#define _SOLVER_RESULTS_H_

#include <string>

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

    void SetSolverName( string& name );
    string& GetSolverName( );

    void SetFitStatisticType( int fitStatType );
    int GetFitStatisticType( );
    
    void StoreBestfitStatisticValue( double fitStatValue );
    double GetBestfitStatisticValue( );

    void StoreNFunctionEvals( int nFunctionEvals );
    int GetNFunctionEvals( );
    
    bool ErrorsPresent( );
    void StoreErrors( double *errors, int nParams );
    void GetErrors( double *errors );


  private:
    int  whichSolver;
    int  whichFitStatistic;
    int  nParameters;
    int  nFuncEvals;
    string  solverName;
    double  bestFitValue;
    bool  paramSigmasPresent;
    bool  paramSigmasAllocated;
    double  *paramSigmas;
    bool  mpResultsPresent;
    mp_result mpResult;
  
};

#endif   // _SOLVER_RESULTS_H_
