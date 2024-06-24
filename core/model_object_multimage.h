/** @file
 * \brief Class declaration for ModelObjectMultImage
 */
/*   Abstract base class interface definition for model_object.cpp [imfit]
 *
 * This is the abstract base class for 1D and 2D "model" objects.
 * 
 */


// CLASS ModelObjectMultImage [derived class]:

#ifndef _MODEL_OBJ_MULTIMAGE_H_
#define _MODEL_OBJ_MULTIMAGE_H_

#include <vector>
#include <string>
#include <tuple>

#include "definitions.h"
#include "function_objects/function_object.h"
#include "model_object.h"
#include "convolver.h"
#include "oversampled_region.h"
#include "param_struct.h"
#include "param_holder.h"
#include "setup_model_object.h"

using namespace std;


// assume that all methods are "common" to both 2D (base) and 1D (derived) versions,
// unless otherwise stated

/// \brief Main class holding data, model information, and code for generating
///        model images, computing chi^2, etc.
class ModelObjectMultImage : public ModelObject
{
  public:
    ModelObjectMultImage( );

    virtual ~ModelObjectMultImage();

    // New methods
    int SetupGeneralOptions( shared_ptr<OptionsBase> options );
    
    void AddModelObject( ModelObject* newModelObj );
    
    std::tuple<int, int> GetImageDimensions( int imageNumber );

    void CreateAllModelImages( double params[] );
    
    void SetupParamHolder( );

    void GetParameterStringsForOneImage( vector<string> &strings, double params[], 
    									int imageNumber, double errors[]=nullptr );

    void GetPerImageFuncStrings( vector<string> &strings, double params[], 
    									int imageNumber, double errors[]=nullptr );

    void GetAllImageDescriptionParameters( double params[], double outputImageParams[] );

    std::tuple<int, int> GetOffsetsForImage( int imageNumber );
    
    string GetDataFilename( int imageNumber );


    // redefined method/member functions (from ModelObject):
    int FinalModelSetup( ) override;
    
    long GetNValidPixels( ) override;
    
    int AddFunction( FunctionObject *newFunctionObj_ptr, bool isGlobalFunc=true ) override;
    
    int FinalSetupForFitting( ) override;
    
    void CreateModelImage( double params[] ) override;

    double GetFitStatistic( double params[] ) override;

    void ComputeDeviates( double yResults[], double params[] ) override;

    // note that since this has a different type signature from the parent
    // class method, we aren't technically overriding
    double * GetModelImageVector( int imageNumber );

    void PrintDescription( ) override;

    int PrintModelParamsToStrings( vector<string> &stringVector, double params[], 
									double errs[], const char *prefix, 
									bool printLimits ) override;

    int GetNImages( ) override;

    int UseBootstrap( ) override;
    
    int MakeBootstrapSample( ) override;

    string PrintModelParamsHorizontalString( const double params[], 
    											const string& separator="\t" ) override;

    string GetParamHeader( ) override;


  private:
  	vector<ModelObject *> modelObjectsVect;
  	int  nModelObjects;
  	int  nImageDescriptionParams;  // sum of image-descr. parameters (for 2nd & subsequent images)
  	int  nGlobalFuncParams;  // sum of global-model function parameters and function-set X0,Y0
  	// NOTE ModelObject's nPerImageParams is used in this class to keep track
  	// of *all* per-image params, from all component ModelObject instances
  	// NOTE: ModelObject's nParamsTot = nImageDescriptionParams + nGlobalFuncParams + nPerImageParams
  	
  	vector<double *> paramsForModelObjects;   // one array-of-double for each ModelObject
  	vector<int> nParamsForModelObjects;   // number of parameters [excluding image-description
  	                 // parameters] for each ModelObject
  	vector<int> nPerImageParamsForModelObjects;  // number of per-image function 
  					 // parameters for each ModelObject [= 0 for those w/o per-image functions]
  	bool paramsForModelObjects_allocated;
  	vector<double *> imageDescriptionParams;
  	bool imageDescriptionParams_allocated;
  	ParamHolder parameterHolder;
  	bool parameterHolderSet;
  
};

#endif   // _MODEL_OBJ_MULTIMAGE_H_
