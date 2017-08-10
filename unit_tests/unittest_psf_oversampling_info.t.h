// Unit tests for PsfOversamplingInfo class
//

// See run_unittest_psf_oversampling_info.sh for how to compile and run these tests.


#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
#include <tuple>
using namespace std;

#include "psf_oversampling_info.h"


double psfPixels0[9] = {0.1, 0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.1, 0.1};
const int nColsPsf = 3;
const int nRowsPsf = 3;
const int nPixelsPsf = 9;

string regionString0 = "10:20,10:20";
string regionString1 = "100:120,180:200";



class TestPsfOversamplingInfo : public CxxTest::TestSuite 
{

public:
  void testFullConstructor( void )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    double *pixVect;
    int outputScale, n;
    int x0_offset, y0_offset;
    string outputString;
    int scale0 = 3;
    
    // allocate using malloc bcs PsfOversamplingInfo destructor will want to call free
    double * psfPixels;
    psfPixels = (double *)malloc(9*sizeof(double));
    for (int i = 0; i < 9; i++)
      psfPixels[i] = psfPixels0[i];
    
    osampleInfo_ptr = new PsfOversamplingInfo(psfPixels, nColsPsf, nRowsPsf, scale0,
    											regionString0);

    n = osampleInfo_ptr->GetNColumns();
    TS_ASSERT_EQUALS(nColsPsf, n);
    n = osampleInfo_ptr->GetNRows();
    TS_ASSERT_EQUALS(nRowsPsf, n);
    pixVect = osampleInfo_ptr->GetPsfPixels();
    for (int i = 0; i < nPixelsPsf; i++)
      TS_ASSERT_EQUALS(pixVect[i], psfPixels0[i]);

    outputScale = osampleInfo_ptr->GetOversamplingScale();
    TS_ASSERT_EQUALS(outputScale, scale0);

    outputString = osampleInfo_ptr->GetRegionString();
    TS_ASSERT_EQUALS(outputString, regionString0);
    
    osampleInfo_ptr->GetImageOffset(x0_offset, y0_offset);
    TS_ASSERT_EQUALS(x0_offset, 0);
    TS_ASSERT_EQUALS(y0_offset, 0);

    delete osampleInfo_ptr;  
  }

  void testFullConstructor_withOffset( void )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    double *pixVect;
    int outputScale, n;
    int x0_offset, y0_offset;
    int x0_offset_in = 100;
    int y0_offset_in = 200;
    string outputString;
    int scale0 = 5;

    // allocate using malloc bcs PsfOversamplingInfo destructor will want to call free
    double * psfPixels;
    psfPixels = (double *)malloc(9*sizeof(double));
    for (int i = 0; i < 9; i++)
      psfPixels[i] = psfPixels0[i];

    osampleInfo_ptr = new PsfOversamplingInfo(psfPixels, nColsPsf, nRowsPsf, scale0,
    											regionString0, x0_offset_in, y0_offset_in);

    n = osampleInfo_ptr->GetNColumns();
    TS_ASSERT_EQUALS(nColsPsf, n);
    n = osampleInfo_ptr->GetNRows();
    TS_ASSERT_EQUALS(nRowsPsf, n);
    pixVect = osampleInfo_ptr->GetPsfPixels();
    for (int i = 0; i < nPixelsPsf; i++)
      TS_ASSERT_EQUALS(pixVect[i], psfPixels0[i]);

    outputScale = osampleInfo_ptr->GetOversamplingScale();
    TS_ASSERT_EQUALS(outputScale, scale0);

    outputString = osampleInfo_ptr->GetRegionString();
    TS_ASSERT_EQUALS(outputString, regionString0);
    
    osampleInfo_ptr->GetImageOffset(x0_offset, y0_offset);
    TS_ASSERT_EQUALS(x0_offset, x0_offset_in);
    TS_ASSERT_EQUALS(y0_offset, y0_offset_in);

    delete osampleInfo_ptr;  
  }
  
  
  void testAddAndRetrieveRegion( void )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    osampleInfo_ptr = new PsfOversamplingInfo();
    string outputString;
    
    
    osampleInfo_ptr->AddRegionString(regionString0);
    outputString = osampleInfo_ptr->GetRegionString();
    TS_ASSERT_EQUALS(outputString, regionString0);

    osampleInfo_ptr->AddRegionString(regionString1);
    outputString = osampleInfo_ptr->GetRegionString();
    TS_ASSERT_EQUALS(outputString, regionString1);
    
    delete osampleInfo_ptr;  
  }

  void testAddAndRetrieveScale( void )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    osampleInfo_ptr = new PsfOversamplingInfo();
    int outputScale;
    int scale0 = 3;
    int scale1 = 5;
    
    
    osampleInfo_ptr->AddOversamplingScale(scale0);
    outputScale = osampleInfo_ptr->GetOversamplingScale();
    TS_ASSERT_EQUALS(outputScale, scale0);

    osampleInfo_ptr->AddOversamplingScale(scale1);
    outputScale = osampleInfo_ptr->GetOversamplingScale();
    TS_ASSERT_EQUALS(outputScale, scale1);
    
    delete osampleInfo_ptr;  
  }

  void testAddAndRetrievePixels( void )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    osampleInfo_ptr = new PsfOversamplingInfo();
    double *pixVect;
    int  n;    
    
    // allocate using malloc bcs PsfOversamplingInfo destructor will want to call free
    double * psfPixels;
    psfPixels = (double *)malloc(9*sizeof(double));
    for (int i = 0; i < 9; i++)
      psfPixels[i] = psfPixels0[i];

    osampleInfo_ptr->AddPsfPixels(psfPixels, nColsPsf, nRowsPsf, true);
    n = osampleInfo_ptr->GetNColumns();
    TS_ASSERT_EQUALS(nColsPsf, n);
    n = osampleInfo_ptr->GetNRows();
    TS_ASSERT_EQUALS(nRowsPsf, n);

    pixVect = osampleInfo_ptr->GetPsfPixels();
    for (int i = 0; i < nPixelsPsf; i++)
      TS_ASSERT_EQUALS(pixVect[i], psfPixels0[i]);

    delete osampleInfo_ptr;  
  }

  void testAddAndRetrieveImageOffset( void )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    osampleInfo_ptr = new PsfOversamplingInfo();
    int outputScale;
    int scale0 = 3;
    int scale1 = 5;
    
    
    osampleInfo_ptr->AddOversamplingScale(scale0);
    outputScale = osampleInfo_ptr->GetOversamplingScale();
    TS_ASSERT_EQUALS(outputScale, scale0);

    osampleInfo_ptr->AddOversamplingScale(scale1);
    outputScale = osampleInfo_ptr->GetOversamplingScale();
    TS_ASSERT_EQUALS(outputScale, scale1);
    
    delete osampleInfo_ptr;  
  }

  void testGetCorrectedRegionCoords1( )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    int x0_offset_in = 0;
    int y0_offset_in = 0;
    int scale0 = 3;
    int  x1, x2, y1, y2;
    int  x1_true, x2_true, y1_true, y2_true;
    string regionString = "101:110,201:210";
    x1_true = 101;
    x2_true = 110;
    y1_true = 201;
    y2_true = 210;

    // allocate using malloc bcs PsfOversamplingInfo destructor will want to call free
    double * psfPixels;
    psfPixels = (double *)malloc(9*sizeof(double));
    for (int i = 0; i < 9; i++)
      psfPixels[i] = psfPixels0[i];

    osampleInfo_ptr = new PsfOversamplingInfo(psfPixels, nColsPsf, nRowsPsf, scale0,
    											regionString, x0_offset_in, y0_offset_in);

    std::tie(x1, x2, y1, y2) = osampleInfo_ptr->GetCorrectedRegionCoords();
    TS_ASSERT_EQUALS(x1, x1_true);
    TS_ASSERT_EQUALS(x2, x2_true);
    TS_ASSERT_EQUALS(y1, y1_true);
    TS_ASSERT_EQUALS(y2, y2_true);

    delete osampleInfo_ptr;  
  }

  void testGetCorrectedRegionCoords2( )
  {
    PsfOversamplingInfo *osampleInfo_ptr;
    int x0_offset_in = 100;
    int y0_offset_in = 200;
    int scale0 = 3;
    int  x1, x2, y1, y2;
    int  x1_true, x2_true, y1_true, y2_true;
    string regionString = "101:110,201:220";
    x1_true = 1;
    x2_true = 10;
    y1_true = 1;
    y2_true = 20;

    // allocate using malloc bcs PsfOversamplingInfo destructor will want to call free
    double * psfPixels;
    psfPixels = (double *)malloc(9*sizeof(double));
    for (int i = 0; i < 9; i++)
      psfPixels[i] = psfPixels0[i];

    osampleInfo_ptr = new PsfOversamplingInfo(psfPixels, nColsPsf, nRowsPsf, scale0,
    											regionString, x0_offset_in, y0_offset_in);

    std::tie(x1, x2, y1, y2) = osampleInfo_ptr->GetCorrectedRegionCoords();
    TS_ASSERT_EQUALS(x1, x1_true);
    TS_ASSERT_EQUALS(x2, x2_true);
    TS_ASSERT_EQUALS(y1, y1_true);
    TS_ASSERT_EQUALS(y2, y2_true);

    delete osampleInfo_ptr;  
  }
};

