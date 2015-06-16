// $CXXTESTGEN --error-printer -o test_runner_solver_results.cpp unit_tests/unittest_solver_results.t.h
// $CPP -std=c++11 -o test_runner_solver_results test_runner_solver_results.cpp solver_results.cpp -I. -I/usr/local/include -I$CXXTEST
// ./test_runner_solver_results

// NOTE: we compile with the -std=c++11 flag in order to make the initialization of
// the internal orig_errors array simpler...

#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
using namespace std;
#include "solver_results.h"
#include "definitions.h"

const double BESTFIT_VAL_MPFIT = 56975.472675;
const int N_FUNCEVALS_MPFIT = 273;


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  SolverResults *solverResults_mpfit;
  SolverResults *solverResults;
  mp_result mpResult_ref;
  // following syntax requires -std=c++11 while compiling
  double orig_errors[16] = {0.0188, 0.0209, 0.839427, 0.0149152, 0.0, 2.1084,
  							0.269065, 0.0429224, 0.00122636, 6.24293, 0.247773,
  							0.0716676, 0.00130667, 1.25109, 0.807103, 0.539445};
  

  void setUp()
  {
    solverResults = new SolverResults;
    solverResults_mpfit = new SolverResults;

    mpResult_ref.bestnorm = BESTFIT_VAL_MPFIT;
    mpResult_ref.orignorm = 267463.538419;
    mpResult_ref.niter = 17;
    mpResult_ref.nfev = N_FUNCEVALS_MPFIT;
    mpResult_ref.status = 1;
    mpResult_ref.npar = 16;
    mpResult_ref.nfree = 16;
    mpResult_ref.npegged = 1;
    mpResult_ref.nfunc = 38213;
    mpResult_ref.resid = NULL;
    mpResult_ref.xerror = orig_errors;
    mpResult_ref.covar = NULL;
  }

  void tearDown()
  {
    delete solverResults;
    delete solverResults_mpfit;
  }


  void testSetAndGetSolverType( void )
  {
    int  output;
    int  correct = MPFIT_SOLVER;

    output = solverResults->GetSolverType();
    TS_ASSERT_EQUALS(correct, output);

    correct = DIFF_EVOLN_SOLVER;
    solverResults->SetSolverType(DIFF_EVOLN_SOLVER);
    output = solverResults->GetSolverType();
    TS_ASSERT_EQUALS(correct, output);
  }

  void testSetAndGetFitStatistic( void )
  {
    int  output;
    int  correct = FITSTAT_CHISQUARE;

    output = solverResults->GetFitStatisticType();
    TS_ASSERT_EQUALS(correct, output);

    correct = FITSTAT_POISSON_MLR;
    solverResults->SetFitStatisticType(FITSTAT_POISSON_MLR);
    output = solverResults->GetFitStatisticType();
    TS_ASSERT_EQUALS(correct, output);
  }

  void testSetAndGetBestFitValue( void )
  {
    double  output;
    double  correct1 = 0.0;
    double  correct2 = 56975.472675;

    output = solverResults->GetBestfitStatisticValue();
    TS_ASSERT_EQUALS(correct1, output);

    solverResults->StoreBestfitStatisticValue(correct2);
    output = solverResults->GetBestfitStatisticValue();
    TS_ASSERT_EQUALS(correct2, output);
  }

  void testGetSetNFunctionEvals( void )
  {
    int  output;
    int  correct1 = 0;    // default value
    int  correct2 = 140;

    output = solverResults->GetNFunctionEvals();
    TS_ASSERT_EQUALS(correct1, output);

    solverResults->StoreNFunctionEvals(correct2);
    output = solverResults->GetNFunctionEvals();
    TS_ASSERT_EQUALS(correct2, output);
  }


  void testAddMPResults( void )
  {
    int output1;
    double output2;
    
    solverResults_mpfit->AddMPResults(mpResult_ref);
 
    output1 = solverResults_mpfit->GetNFunctionEvals();
    TS_ASSERT_EQUALS(N_FUNCEVALS_MPFIT, output1);

    output2 = solverResults_mpfit->GetBestfitStatisticValue();
    TS_ASSERT_EQUALS(BESTFIT_VAL_MPFIT, output2);
   
  }

};
