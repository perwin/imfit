/*   Abstract base lass interface definition for model_object.cpp [imfit]
 *   VERSION 0.3
 *
 * This is the abstract base clase for 1D and 2D "model" objects.
 * 
 */


// CLASS ModelObject [base class]:

#ifndef _MODEL_OBJ_H_
#define _MODEL_OBJ_H_

#include <vector>

#include "definitions.h"
#include "function_object.h"
#include "convolver.h"
#include "param_struct.h"

using namespace std;


class ModelObject
{
  public:
    // Constructors:
    ModelObject( );

    void SetDebugLevel( int debuggingLevel );
    
    // common, not specialized
    void AddFunction( FunctionObject *newFunctionObj_ptr );
    
    // common, but Specialized by ModelObject1D
    virtual void DefineFunctionSets( vector<int>& functionStartIndices );
    
    // 1D only, but needs to be part of base interface
    virtual void AddDataVectors( int nDataValues, double *xValVector, 
    						double *yValVector, bool magnitudeData );

		// 2D only
    void AddImageDataVector( double *pixelVector, int nImageColumns, int nImageRows,
    													int nCombinedImages );

		// 2D only
    void SetupModelImage( int nDataValues, int nImageColumns, int nImageRows );
    
		// 2D only
    virtual void AddErrorVector( int nDataValues, int nImageColumns, int nImageRows,
                         double *pixelVector, int inputType );

    // 1D only
    virtual void AddErrorVector1D( int nDataValues, double *pixelVector, int inputType );

    // 1D only
    virtual void AddMaskVector1D( int nDataValues, double *inputVector, int inputType );

		// 2D only
    virtual void GenerateErrorVector( double gain, double readNoise, double skyValue );

		// 2D only
    virtual void AddMaskVector( int nDataValues, int nImageColumns, int nImageRows,
                         double *pixelVector, int inputType );

		// 2D only
    void AddPSFVector( int nPixels_psf, int nColumns_psf, int nRows_psf,
                         double *psfPixels );

    // 1D only
    virtual void AddPSFVector1D( int nPixels_psf, double *xValVector, double *yValVector );
    
		// 2D only [1D maybe needs something similar, but with diff. interface]
    virtual void ApplyMask( );

    // common, but Specialized by ModelObject1D
    virtual void CreateModelImage( double params[] );
    
    // Specialized by ModelObject1D
    virtual void ComputeDeviates( double yResults[], double params[] );

    // common, not specialized
    virtual void SetupChisquaredCalcs( );
    
    // common, not specialized
    virtual double ChiSquared( double params[] );
    
    // common, but Specialized by ModelObject1D
    virtual void PrintDescription( );

    // common, but Specialized by ModelObject1D
    virtual void PrintModelParams( FILE *output_ptr, double params[], mp_par *parameterInfo,
																		double errs[] );

    // 2D only; NOT USED ANYWHERE!
    void PrintImage( double *pixelVector );

		// 2D only
    void PrintInputImage( );

		// 2D only
    void PrintModelImage( );

    // 2D only; NOT USED ANYWHERE!
    void PrintWeights( );

    // common, but Specialized by ModelObject1D
    virtual void PopulateParameterNames( );

    // common, might be specialized...
    virtual void FinalSetup( );

    // common, not specialized
    string& GetParameterName( int i );

    // common, not specialized
    int GetNParams( );

    // common, not specialized
    int GetNValidPixels( );

		// 2D only
    double * GetModelImageVector( );

    // 1D only
    virtual int GetModelVector( double *profileVector );

    // Destructor
    virtual ~ModelObject();


  private:
    Convolver  *psfConvolver;
  
  protected:  // same as private, except accessible to derived classes
    int  nDataVals, nColumns, nRows, nValidDataVals, nCombined;
    double  nCombined_sqrt;
    int  debugLevel;
    bool  dataValsSet, parameterBoundsSet, modelVectorAllocated, weightVectorAllocated;
    bool  setStartFlag_allocated;
    bool  modelImageComputed;
    bool  weightValsSet, maskExists;
    bool  doConvolution;
    bool  doChisquared;   // ModelObject will be asked to do chi-squared calculations
    int  nFunctions, nFunctionSets, nFunctionParams, nParamsTot;
    double  *dataVector;
    double  *weightVector;
    double  *maskVector;
    double  *modelVector;
    double  *deviatesVector;
    double  *parameterBounds;
    int  *functionSetStarts;
    bool  *setStartFlag;
    vector<FunctionObject *> functionObjects;
    vector<int> paramSizes;
    vector<string>  parameterLabels;
    
    bool CheckParamVector( int nParams, double paramVector[] );
    bool CheckWeightVector( );
    bool VetDataVector( );
  
};

#endif   // _MODEL_OBJ_H_
