#ifndef _PSF_INTERPOLATORS_H_
#define _PSF_INTERPOLATORS_H_

#include "gsl/gsl_spline2d.h"

#define kInterpolator_Base 0
#define kInterpolator_bicubic 1
#define kInterpolator_lanczos2 2
#define kInterpolator_lanczos3 3


class PsfInterpolator
{
  public:
  // need to provide zero-parameter base-class constructor, since derived-class
  // constructors will automatically try to call it
  PsfInterpolator( ) { ; };
  // derived classes should implement their own versions of this
  PsfInterpolator( double *inputImage, int nCols_image, int nRows_image ) { ; };
  
  virtual ~PsfInterpolator( ) { ; };
  
  int GetInterpolatorType( ) { return interpolatorType; };
  
  virtual double GetValue( double x, double y ) = 0;

  protected:
    // class constant
    const static int  interpolatorType = kInterpolator_Base;
    // data members proper
    int  nColumns, nRows;
    long  nPixelsTot;
    double  xBound, yBound, deltaXMin, deltaXMax, deltaYMin, deltaYMax;
};


// Derived class using GNU Scientific Library's 2D bicubic interpolation
class PsfInterpolator_bicubic : public PsfInterpolator
{
  public:
  PsfInterpolator_bicubic( double *inputImage, int nCols_image, int nRows_image );
  
  ~PsfInterpolator_bicubic( );
  
  double GetValue( double x, double y );

  protected:
    // class constant
    const static int  interpolatorType = kInterpolator_bicubic;
    
  private:
    // new data members
    gsl_spline2d *splineInterp;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    double *xArray;
    double *yArray;
};

#endif   // _PSF_INTERPOLATORS_H_
