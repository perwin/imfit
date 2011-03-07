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


// CLASS ModelObject1d [base class]:

#ifndef _MODEL_OBJ_1D_H_
#define _MODEL_OBJ_1D_H_

#include <vector>

#include "definitions.h"
#include "function_object.h"
#include "model_object.h"
#include "convolver1d.h"

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

    void AddMaskVector1D( int nDataValues, double *inputVector, int inputType );
    
    void AddPSFVector1D( int nPixels_psf, double *xValVector, double *yValVector );
    
    void CreateModelImage( double params[] );

    void ComputeDeviates( double yResults[], double params[] );

    void PrintDescription( );
    
    void PrintModelParams( FILE *output_ptr, double params[], mp_par *parameterInfo,
																		double errs[] );
    void PopulateParameterNames( );
    
    int GetModelVector( double *profileVector );

    void UseBootstrap( );
    
    void MakeBootstrapSample( );

    // Destructor
    ~ModelObject1d();


  private:
    Convolver1D  *psfConvolver;
  
//   protected:  // same as private, except accessible to derived classes
//     int  nParams;
//     int  nDataVals, nColumns, nRows;
//     bool  dataValsSet, parameterBoundsSet, modelVectorAllocated, modelImageComputed;
//     bool  weightValsSet;
//     int  nFunctions, nParamsTot;
// 	double  *dataVector;
	  double  *dataXValues, *modelXValues;
	  double  zeroPoint;
	  bool  dataAreMagnitudes, doBootstrap, bootstrapIndicesAllocated;
	// things useful for PSF convolution
	  int  dataStartOffset, nPSFVals, nModelVals;
	  int  *bootstrapIndices;
// 	double  *weightVector;
// 	double  *modelVector;
// 	double  *parameterBounds;
// 	std::vector<FunctionObject *> functionObjects;
// 	std::vector<int> paramSizes;
  
};

#endif   // _MODEL_OBJ_1D_H_
