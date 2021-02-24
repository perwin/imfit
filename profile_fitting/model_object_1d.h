/*   Class interface definition for model_object_1d.cpp
 *   VERSION 0.1
 *
 * This is intended to be an abstract base class for the various
 * "model" objects (e.g., image data + fitting functions).
 * 
 *   We wrap the declaration in an ifndef -- endif pair to
 * prevent this from being read in more than once (e.g., from
 * *main.cpp files which read in header files for more than
 * one derived class).
 */


// CLASS ModelObject1d [derived class]:

#ifndef _MODEL_OBJ_1D_H_
#define _MODEL_OBJ_1D_H_

#include <vector>

#include "definitions.h"
#include "function_object.h"
#include "model_object.h"
#include "convolver1d.h"
#include "param_struct.h"

class ModelObject1d : public ModelObject
{
  public:
    // Constructors:
    ModelObject1d( );
    
   // redefined method/member functions:
    void DefineFunctionSets( vector<int>& functionStartIndices );
    
    void AddDataVectors( int nDataValues, double *xValVector, double *yValVector,
    											bool magnitudeData );

    void SetZeroPoint( double zeroPointValue );
    
    void AddErrorVector1D( int nDataValues, double *inputVector, int inputType );

    int AddMaskVector1D( int nDataValues, double *inputVector, int inputType );
    
    int AddPSFVector1D( int nPixels_psf, double *xValVector, double *yValVector );
    
    void CreateModelImage( double params[] );

    void ComputeDeviates( double yResults[], double params[] );

    void PrintDescription( );
    
    int Dimensionality( ) { return 1;};
    
//     void PrintModelParams( FILE *output_ptr, double params[], double errs[], 
//     						const char *prefix="" );

    int PrintModelParamsToStrings( vector<string> &stringVector, double params[], 
									double errs[], const char *prefix="", 
									bool printLimits=false ) override;

    void PopulateParameterNames( );
    
    int FinalSetupForFitting( );

    int GetModelVector( double *profileVector );

//     int UseBootstrap( );
    
//     int MakeBootstrapSample( );

    void PrintVector( double *theVector, int nVals );
    
    void PrintInputImage( );

    void PrintModelImage( );

    void PrintMask( );

    void PrintWeights( );

    // Destructor
    ~ModelObject1d();


  private:
    Convolver1D  *psfConvolver;
  
      bool  parameterBoundsSet;
	  double  *dataXValues, *modelXValues;
	  bool  dataAreMagnitudes;
	// things useful for PSF convolution
	  int  dataStartOffset, nPSFVals;
	  double  *parameterBounds;
};

#endif   // _MODEL_OBJ_1D_H_
