/*   Class interface definition for model_object.cpp [imfit]
 *   VERSION 0.3
 *
 * This is intended to be an abstract base class for the various
 * "model" objects (e.g., image data + fitting functions).
 * 
 */


// CLASS ModelObject [base class]:

#ifndef _MODEL_OBJ_H_
#define _MODEL_OBJ_H_

#include <vector>

#include "definitions.h"
#include "function_object.h"

using namespace std;


class ModelObject
{
  public:
    // Constructors:
    ModelObject( );
    
    void AddFunction( FunctionObject *newFunctionObj_ptr );
    
    void DefineFunctionSets( vector<int>& functionStartIndices );
    
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

    virtual void ApplyMask( );
    
    void AddParameterLimits( double *paramLimits );
    
    virtual void CreateModelImage( double params[] );

    // may need to be overridden by derived class ModelObject1D
    virtual void ComputeDeviates( double yResults[], double params[] );

    virtual double ChiSquared( double params[] );
    
    virtual void GetDescription( );

    void PrintImage( double *pixelVector );

    void PrintInputImage( );

    void PrintModelImage( );

    void PrintWeights( );

    void PopulateParameterNames( );

    string& GetParameterName( int i );

    int GetNParams( );

    int GetNValidPixels( );

    double * GetModelImageVector( );

    // Destructor
    virtual ~ModelObject();


  private:
  
  protected:  // same as private, except accessible to derived classes
    int  nDataVals, nColumns, nRows, nValidDataVals;
    bool  dataValsSet, parameterBoundsSet, modelVectorAllocated, weightVectorAllocated;
    bool  modelImageComputed;
    bool  weightValsSet, maskExists;
    int  nFunctions, nFunctionSets, nFunctionParams, nParamsTot;
    double  *dataVector;
    double  *weightVector;
    double  *maskVector;
    double  *modelVector;
    double  *parameterBounds;
    int  *functionSetStarts;
    bool  *setStartFlag;
    vector<FunctionObject *> functionObjects;
    vector<int> paramSizes;
    vector<string>  parameterLabels;
    
    bool CheckWeightVector( );
  
};

#endif   // _MODEL_OBJ_H_
