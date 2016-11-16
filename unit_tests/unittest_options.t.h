// Unit tests for options classes
//

// See run_unittest_options.sh for how to compile and run these tests.


#include <cxxtest/TestSuite.h>

#include <math.h>
#include <string>
#include <vector>
using namespace std;

#include "definitions.h"
#include "options_base.h"
#include "options_makeimage.h"
#include "options_imfit.h"
#include "options_mcmc.h"

class TestMakeimageOptions : public CxxTest::TestSuite 
{

public:
  void testBasic( void )
  {
    OptionsBase *baseOptions_ptr;
    
    baseOptions_ptr = new OptionsBase();
    
//    TS_ASSERT_EQUALS( baseOptions_ptr->noConfigFile, true );
    TS_ASSERT_EQUALS( baseOptions_ptr->configFileName, string("") );

    delete baseOptions_ptr;
  }

  void testMakeimageOptions( void )
  {
    OptionsBase *baseOptions_ptr;
    MakeimageOptions *makeimageOptions_ptr;
    
    baseOptions_ptr = new MakeimageOptions();
    
    TS_ASSERT_EQUALS( baseOptions_ptr->noConfigFile, true );
    TS_ASSERT_EQUALS( baseOptions_ptr->configFileName, string("") );

    // pointer to MakeimageOptions object that compiler will recognize as such
    makeimageOptions_ptr = (MakeimageOptions *)baseOptions_ptr;
    TS_ASSERT_EQUALS( makeimageOptions_ptr->noOutputImageName, true );
    TS_ASSERT_EQUALS( makeimageOptions_ptr->configFileName, "" );
    TS_ASSERT_EQUALS( makeimageOptions_ptr->noOutputImageName, true );
    TS_ASSERT_EQUALS( makeimageOptions_ptr->outputImageName, DEFAULT_MAKEIMAGE_OUTPUT_FILENAME );
    TS_ASSERT_EQUALS( makeimageOptions_ptr->noRefImage, true );
    TS_ASSERT_EQUALS( makeimageOptions_ptr->referenceImageName, "" );
    
    TS_ASSERT_EQUALS( makeimageOptions_ptr->printFluxes, false );
    TS_ASSERT_EQUALS( makeimageOptions_ptr->estimationImageSize, DEFAULT_ESTIMATION_IMAGE_SIZE );
    TS_ASSERT_EQUALS( makeimageOptions_ptr->timingIterations, 0 );

    delete baseOptions_ptr;
  }

  void testImfitOptions( void )
  {
    OptionsBase *baseOptions_ptr;
    ImfitOptions *imfitOptions_ptr;
    
    baseOptions_ptr = new ImfitOptions();
    
    TS_ASSERT_EQUALS( baseOptions_ptr->noConfigFile, true );
    TS_ASSERT_EQUALS( baseOptions_ptr->configFileName, string(DEFAULT_IMFIT_CONFIG_FILE) );
    TS_ASSERT_EQUALS( baseOptions_ptr->gain, 1.0 );

    // pointer to ImfitOptions object that compiler will recognize as such
    imfitOptions_ptr = (ImfitOptions *)baseOptions_ptr;
    TS_ASSERT_EQUALS( imfitOptions_ptr->outputParameterFileName, DEFAULT_OUTPUT_PARAMETER_FILE );
    TS_ASSERT_EQUALS( imfitOptions_ptr->gainSet, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->gain, 1.0 );
    TS_ASSERT_EQUALS( imfitOptions_ptr->readNoiseSet, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->readNoise, 0.0 );
    TS_ASSERT_EQUALS( imfitOptions_ptr->expTimeSet, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->expTime, 1.0 );
    TS_ASSERT_EQUALS( imfitOptions_ptr->nCombinedSet, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->nCombined, 1 );
    TS_ASSERT_EQUALS( imfitOptions_ptr->originalSkySet, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->originalSky, 0.0 );
    
    TS_ASSERT_EQUALS( imfitOptions_ptr->useModelForErrors, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->useCashStatistic, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->usePoissonMLR, false );

    TS_ASSERT_EQUALS( imfitOptions_ptr->noParamLimits, true );
    TS_ASSERT_EQUALS( imfitOptions_ptr->subsamplingFlag, true );

    TS_ASSERT_EQUALS( imfitOptions_ptr->doBootstrap, false );
    TS_ASSERT_EQUALS( imfitOptions_ptr->bootstrapIterations, 0 );
    TS_ASSERT_EQUALS( imfitOptions_ptr->rngSeed, 0 );

    delete baseOptions_ptr;
  }

  void testMCMCOptions( void )
  {
    OptionsBase *baseOptions_ptr;
    MCMCOptions *mcmcOptions_ptr;
    
    baseOptions_ptr = new MCMCOptions();
    
    TS_ASSERT_EQUALS( baseOptions_ptr->noConfigFile, true );
    TS_ASSERT_EQUALS( baseOptions_ptr->configFileName, string(DEFAULT_IMFIT_CONFIG_FILE) );
    TS_ASSERT_EQUALS( baseOptions_ptr->gain, 1.0 );

    // pointer to MCMCOptions object that compiler will recognize as such
    mcmcOptions_ptr = (MCMCOptions *)baseOptions_ptr;
    TS_ASSERT_EQUALS( mcmcOptions_ptr->gainSet, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->gain, 1.0 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->readNoiseSet, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->readNoise, 0.0 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->expTimeSet, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->expTime, 1.0 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->nCombinedSet, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->nCombined, 1 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->originalSkySet, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->originalSky, 0.0 );
    
    TS_ASSERT_EQUALS( mcmcOptions_ptr->useModelForErrors, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->useCashStatistic, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->usePoissonMLR, false );

    TS_ASSERT_EQUALS( mcmcOptions_ptr->noParamLimits, true );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->subsamplingFlag, true );

    TS_ASSERT_EQUALS( mcmcOptions_ptr->appendToOutput, false );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->outputFileRoot, "mcmc_out" );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->nChains, -1 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->maxEvals, 100000 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->nBurnIn, 5000 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->GRScaleReductionLimit, 1.01 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->mcmcNoise, 0.01 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->mcmc_bstar, 1.0e-6 );
    TS_ASSERT_EQUALS( mcmcOptions_ptr->rngSeed, 0 );

    delete baseOptions_ptr;
  }

};

