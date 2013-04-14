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
#include <string>

#include "definitions.h"
#include "function_objects/function_object.h"
#include "convolver.h"
#include "param_struct.h"

using namespace std;


class ModelObject
{
  public:
    // Constructors:
    ModelObject( );

    void SetDebugLevel( int debuggingLevel );
    
    void SetMaxThreads( int maxThreadNumber );
    
    // common, not specialized
    // Adds a new FunctionObject pointer to the internal vector
    void AddFunction( FunctionObject *newFunctionObj_ptr );
    
    // common, but Specialized by ModelObject1D
    virtual void DefineFunctionSets( vector<int>& functionStartIndices );
    
    // 1D only, but needs to be part of base interface
    virtual void AddDataVectors( int nDataValues, double *xValVector, 
    						double *yValVector, bool magnitudeData ) { nDataVals = nDataValues; };

    // Probably 1D only, but might be usable by 2D version later...
    virtual void SetZeroPoint( double zeroPointValue );

 
//     void SetGain( double gainValue );
// 
//     void SetSkyBackground( double originalSkyValue );

	// 2D only
    void AddImageDataVector( double *pixelVector, int nImageColumns, int nImageRows );

	// 2D only
    void AddImageCharacteristics( double imageGain, double readoutNoise, double expTime, 
    							int nCombinedImages, double originalSkyBackground );
    
	// 2D only
    void SetupModelImage( int nImageColumns, int nImageRows );
    
	// 2D only
    virtual void AddErrorVector( int nDataValues, int nImageColumns, int nImageRows,
                         double *pixelVector, int inputType );

    // 1D only
    virtual void AddErrorVector1D( int nDataValues, double *pixelVector, int inputType ) { ; };

    // 1D only
    virtual void AddMaskVector1D( int nDataValues, double *inputVector, int inputType ) { ; };
    
	// 2D only
    virtual void GenerateErrorVector( );

	// 2D only
    virtual void AddMaskVector( int nDataValues, int nImageColumns, int nImageRows,
                         double *pixelVector, int inputType );

	// 2D only
    void AddPSFVector( int nPixels_psf, int nColumns_psf, int nRows_psf,
                         double *psfPixels );

    // 1D only
    virtual void AddPSFVector1D( int nPixels_psf, double *xValVector, double *yValVector ) { ; };
    
	// 2D only [1D maybe needs something similar, but with diff. interface]
    virtual void ApplyMask( );

    // common, but Specialized by ModelObject1D
    virtual void CreateModelImage( double params[] );
    
    // 2D only
    void UpdateWeightVector(  );

    // Specialized by ModelObject1D
    virtual void ComputeDeviates( double yResults[], double params[] );

     // common, not specialized
    virtual void UseModelErrors( );

     // common, not specialized
    virtual void UseCashStatistic( );
 
     // common, not specialized
    virtual bool UsingCashStatistic( );
 
    // common, not specialized
    virtual double GetFitStatistic( double params[] );
    
    // common, not specialized
    virtual double ChiSquared( double params[] );
    
    // common, not specialized
    virtual double CashStatistic( double params[] );
    
    // common, but Specialized by ModelObject1D
    virtual void PrintDescription( );

    // common, not specialized
    void GetFunctionNames( vector<string>& functionNames );

    // common, but Specialized by ModelObject1D
    virtual void PrintModelParams( FILE *output_ptr, double params[], mp_par *parameterInfo,
																		double errs[] );

    // 2D only; NOT USED ANYWHERE!
    void PrintImage( double *pixelVector, int nColumns, int nRows );

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
    int GetNFunctions( );

    // common, not specialized
    int GetNParams( );

    // common, not specialized -- returns total number of data values (e.g., pixels in image)
    int GetNDataValues( );

    // common, not specialized -- returns total number of *non-masked* data values
    int GetNValidPixels( );

		// 2D only
    double * GetModelImageVector( );

		// 2D only
    double * GetExpandedModelImageVector( );

		// 2D only
    double * GetResidualImageVector( );

		// 2D only
    double * GetWeightImageVector( );

		// 2D only
    double FindTotalFluxes(double params[], int xSize, int ySize, 
    											double individualFluxes[] );

    // Generate a model image using *one* of the FunctionObjects (the one indicated by
    // functionIndex) and the input parameter vector; returns pointer to modelVector.
    double * GetSingleFunctionImage( double params[], int functionIndex );

    // 1D only
    virtual int GetModelVector( double *profileVector ) { return -1; };

    // 1D or 2D
    virtual void UseBootstrap( );
    
    // 1D or 2D
    virtual void MakeBootstrapSample( );
    
    // Destructor
    virtual ~ModelObject();


  private:
    Convolver  *psfConvolver;
  
  protected:  // same as private, except accessible to derived classes
    int  nDataVals, nDataColumns, nDataRows, nValidDataVals, nCombined;
    int  nModelVals, nModelColumns, nModelRows, nPSFColumns, nPSFRows;
//    double  nCombined_sqrt;
	double  zeroPoint;
	double  gain, readNoise, exposureTime, originalSky, effectiveGain;
	double  readNoise_adu_squared;
    int  debugLevel;
    int  maxRequestedThreads;
    bool  dataValsSet, parameterBoundsSet, modelVectorAllocated, weightVectorAllocated;
    bool  residualVectorAllocated, outputModelVectorAllocated;
    bool  setStartFlag_allocated;
    bool  modelImageComputed;
    bool  weightValsSet, maskExists, doBootstrap, bootstrapIndicesAllocated;
    bool  doConvolution;
    bool  modelErrors;
    bool  useCashStatistic;
    bool  deviatesVectorAllocated;   // for chi-squared calculations
    bool  zeroPointSet;
    int  nFunctions, nFunctionSets, nFunctionParams, nParamsTot;
    double  *dataVector;
    double  *weightVector;
    double  *maskVector;
    double  *modelVector;
    double  *deviatesVector;
    double  *residualVector;
    double  *outputModelVector;
    double  *parameterBounds;
    int  *bootstrapIndices;
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
