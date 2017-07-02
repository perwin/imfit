// Header file for PsfOversamplingInfo class

#ifndef _PSF_OVERSAMPLING_INFO_H_
#define _PSF_OVERSAMPLING_INFO_H_

#include <string>
#include <vector>

using namespace std;


class PsfOversamplingInfo
{
  public:
    PsfOversamplingInfo();
    PsfOversamplingInfo( double *inputPixels, int nCols, int nRows, int scale,
    					string inputRegionString );
    ~PsfOversamplingInfo();
  
    void AddPsfPixels( double *inputPixels, int nCols, int nRows );
    void AddRegionString( string inputRegionString );
    void AddOversamplingScale( int scale );
  
    int GetNColumns( );
    int GetNRows( );
    double * GetPsfPixels( );
    string GetRegionString( );
    int GetOversamplingScale( );

  private:
    int  nColumns_psf, nRows_psf;
    string  regionString;
    double *  psfPixels;
    int  oversamplingScale;
};


#endif   // _PSF_OVERSAMPLING_INFO_H_
