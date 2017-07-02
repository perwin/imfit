// Unit tests for PsfOversamplingInfo class
//

// See run_unittest_psf_oversampling_info.sh for how to compile and run these tests.


#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
using namespace std;

#include "psf_oversampling_info.h"


double psfPixels0[9] = {0.1, 0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.1, 0.1};
double psfPixels1[9] = {0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1};
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
    string outputString;
    int scale0 = 3;

    osampleInfo_ptr = new PsfOversamplingInfo(psfPixels0, nColsPsf, nRowsPsf, scale0,
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
    
    osampleInfo_ptr->AddPsfPixels(psfPixels0, nColsPsf, nRowsPsf);
    n = osampleInfo_ptr->GetNColumns();
    TS_ASSERT_EQUALS(nColsPsf, n);
    n = osampleInfo_ptr->GetNRows();
    TS_ASSERT_EQUALS(nRowsPsf, n);

    pixVect = osampleInfo_ptr->GetPsfPixels();
    for (int i = 0; i < nPixelsPsf; i++)
      TS_ASSERT_EQUALS(pixVect[i], psfPixels0[i]);

    delete osampleInfo_ptr;  
  }

};

