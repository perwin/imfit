// Header file for ParamHolder class

#ifndef _PARAM_HOLDER_H_
#define _PARAM_HOLDER_H_

#include <vector>

#include "definitions_multimage.h"
#include "paramvector_processing.h"


class ParamHolder
{
  public:
    ParamHolder( );

    ~ParamHolder();

    void GeneralSetup( int nImagesTot, int nFuncs, int nParamsTot, int nFuncSets,
						const vector<int> paramSizes, bool fsetStartFlags[],
						const vector<int> localParamsCount );
    
    void AddNewParameterVector( double params[] );
    
    /// Returns the current set of image-description parameters for image number nImage,
    /// storing them in imageParamsArray as {pixelScale, imageRotation, intensityScale}
    void GetImageParams( int nImage, double *imageParamsArray );

    /// Returns function parameters with X0,Y0 values adjusted for image number nImage, 
    /// storing them in paramsForModelObject
    void GetParamsForModelObject( int nImage, double *paramsForModelObject );
    
    
    // PRIVATE (declared public so we can do unit tests on it)
    void AssembleParametersForImage( double externalInputParamsVect[], int imageNumber, 
								double outputModelParams[] );


  private:
    bool  arraysAllocated;
    double  *pixelScales;
    double  *rotationAngles;
    double  *intensityScales;
    double  *imageX0s;
    double  *imageY0s;
    double  *modelParamArray;
    bool  perImageModelParamsVect_allocated;
    vector<double *>  perImageModelParameters;
    int  nParametersTot;    // *All* parameters (global + local functions + image-desc.)
    int  nImageParameters;  // = N_IMAGE_PARAMS*(nImagesTot - 1)
    int  nModelParameters;  // = nInputParamsTot - nImageParams
    int  nFunctionSets;
    int  nLocalParams;      // all local (per-image-function) parameters
    int  nImages;
    int  nFunctions;
    bool *functionSetStartFlags;  // whether or not a function is at start of new
    						 		// function block; len = nFunctions
    vector<int>  parameterSizes;  // number of parameters per function; len = nFunctions
    vector<int>  nLocalParamsForImages;  // how many per-image function params for each image
    vector<int>  nParamsPerModelObj;  // # params for each image's ModelObject
  
};

#endif   // _PARAM_HOLDER_H_
