#ifndef _PSF_INTERPOLATOR_H_
#define _PSF_INTERPOLATOR_H_

#include "gsl/gsl_spline2d.h"

class PsfInterpolator
{
  public:
  PsfInterpolator( double *inputImage, int nCols_image, int nRows_image );
  
  ~PsfInterpolator( );
  
  double GetValue( double x, double y );

  private:
    int  nColumns, nRows;
    long  nPixelsTot;
    double  xBound, yBound, deltaXMin, deltaXMax, deltaYMin, deltaYMax;
    gsl_spline2d *splineInterp;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    double *xArray;
    double *yArray;
};

#endif   // _PSF_INTERPOLATOR_H_
