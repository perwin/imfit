/*   Class interface definition for model_object.cpp [imfit]
 *   VERSION 0.3
 *
 * This is the standard class for 2D image-generation and image-fitting
 * "model" objects; it also acts as a base class for variant model objects
 * (e.g., 1-D fitting objects -- see model_object_1d.h)
 * 
 */


// CLASS ModelObject [base class]:

#ifndef _MODEL_OBJ_H_
#define _MODEL_OBJ_H_

#include <vector>

#include "definitions.h"
#include "function_object.h"
#include "convolver.h"

using namespace std;


class ModelObject
{
  public:
    // Constructors:
    ModelObject( );
    
    void AddFunction( FunctionObject *newFunctionObj_ptr );
    
    // may need to be overridden by derived class ModelObject1D
    virtual void DefineFunctionSets( vector<int>& functionStartIndices );
    
    // to be overridden by derived class ModelObject1D
    virtual void AddDataVectors( int nDataValues, double *xValVector, 
    						double *yValVector, bool magnitudeData );

    void AddImageDataVector( int nDataValues, int nImageColumns,
                                     int nImageRows, double *pixelVector );

    void SetupModelImage( int nDataValues, int nImageColumns, int nImageRows );
    
    virtual void AddErrorVector( int nDataValues, int nImageColumns, int nImageRows,
                         double *pixelVector, int inputType );

    virtual void AddErrorVector1D( int nDataValues, double *pixelVector, int inputType );

    virtual void GenerateErrorVector( double gain, double readNoise, double skyValue );

    virtual void AddMaskVector( int nDataValues, int nImageColumns, int nImageRows,
                         double *pixelVector, int inputType );

    void AddPSFVector(int nPixels_psf, int nColumns_psf, int nRows_psf,
                         double *psfPixels);

    virtual void ApplyMask( );
    
    void AddParameterLimits( double *paramLimits );
    
    // may need to be overridden by derived class ModelObject1D
    virtual void CreateModelImage( double params[] );
    
    // may need to be overridden by derived class ModelObject1D
    virtual void ComputeDeviates( double yResults[], double params[] );

    virtual void SetupChisquaredCalcs( );
    
    virtual double ChiSquared( double params[] );
    
    virtual void PrintDescription( );

    void PrintImage( double *pixelVector );

    void PrintInputImage( );

    void PrintModelImage( );

    void PrintWeights( );

    // may need to be overridden by derived class ModelObject1D
    virtual void PopulateParameterNames( );

    string& GetParameterName( int i );

    int GetNParams( );

    int GetNValidPixels( );

    double * GetModelImageVector( );

    // Destructor
    virtual ~ModelObject();


  private:
    Convolver  *psfConvolver;
  
  protected:  // same as private, except accessible to derived classes
    int  nDataVals, nColumns, nRows, nValidDataVals;
    bool  dataValsSet, parameterBoundsSet, modelVectorAllocated, weightVectorAllocated;
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
  
};

#endif   // _MODEL_OBJ_H_
