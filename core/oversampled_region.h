/*! \file
   \brief  Class declaration for OversampledRegion (deals with computing 
           oversampled regions of an image and (optionally) convolving them 
           with an oversampled PSF).
 */
 


#ifndef _OVERSAMPLED_REGION_H_
#define _OVERSAMPLED_REGION_H_

#include <string>
#include <vector>

#include "convolver.h"
#include "function_objects/function_object.h"

using namespace std;




/// \brief Class for computing oversampled model image region and downsampling to match main
///        image, with optional PSF convolution using oversampled PSF
class OversampledRegion
{
  public:
    // Constructors and Destructors:
    OversampledRegion( );
    ~OversampledRegion( );
    
    // Public member functions:
    void SetDebugImageName( const string imageName );
    
    void AddPSFVector( double *psfPixels_input, const int nColumns, const int nRows );
    
    void SetMaxThreads( const int maximumThreadNumber );

    void SetDebugLevel( const int debuggingLevel );

    int SetupModelImage( int x1, int y1, int nBaseColumns, int nBaseRows, 
    					int nColumnsMain, int nRowsMain, int nColumnsPSF_main,
    					int nRowsPSF_main, int oversampScale );
    					
    void ComputeRegionAndDownsample( double *mainImageVector, 
    				vector<FunctionObject *> functionObjectVect, int nFunctionObjects );


  private:
  // Private member functions:

  // Data members:
    Convolver  *psfConvolver;
    int  ompChunkSize, maxRequestedThreads, debugLevel;
    int  oversamplingScale;
    double  subpixFrac, startX_offset, startY_offset;
    int  nPSFColumns, nPSFRows;
    int  nRegionColumns, nRegionRows, nRegionVals;
    int  x1_region, y1_region;
    int  nMainImageColumns, nMainImageRows, nMainPSFColumns, nMainPSFRows;
    int  nModelColumns, nModelRows, nModelVals;
    bool  doConvolution, setupComplete, modelVectorAllocated;
    double  *modelVector;
    string  debugImageName;

};


#endif  // _OVERSAMPLED_REGION_H_
