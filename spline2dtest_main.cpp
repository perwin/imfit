// Test code for 2D spline interpolation


// Things to deal with:
//
//   [X] Make sure that GSL 2D spline interpolation doesn't switch x,y
// 
// From GSL manual:
// with i = 0,...,xsize-1 and j = 0,...,ysize-1
// z_ij = za[j*xsize + i]
// So:
// z[x,y] = z[j*nColumns + i]
// 
// 
// In ModelObject: x is indexed by j, y is indexed by i
// j = k % nModelColumns;
// i = k / nModelColumns;
// y = (double)(i + 1);
// x = (double)(j + 1);
// modelVector[i*nModelColumns + j] = newValSum;

//   [] Internal coordinate system for PSF array vs pixel coordinates
//   -- image function will get X0,Y0 and then x,y
//      positions relative to PSF center: xDiff = x - X0, yDiff = y - Y0
//     
//    Ideal "native" coordinate system has PSF center at (0,0), with coordinates
//    running from x = -nCols/2,..,nCols/2
//      PSF with even number of columns, e.g.., nCols = 6
//        -2.5,-1.5,-0.5,0.5,1.5,2.5
//      PSF with odd number of columns, e.g.., nCols = 5
//        -2,-1,0,1,2
//
//
//    CODE TO GENERATE properly centered xArray and yArray
//      xBound = (nColumns - 1) / 2.0;
//      double *xArray = (double *)calloc(nColumns, sizeof(double));
//      for (int n = 0; n < nColumns; n++) {
//        xArray[n] = n - xBound;
//      yBound = (nRows - 1) / 2.0;
//      double *yArray = (double *)calloc(nRows, sizeof(double));
//      for (int n = 0; n < nRows; n++) {
//        yArray[n] = n - yBound;
// 
//    [] How to handle points beyond boundary of PSF
//      for |xDiff| > halfWidth_psf --> return 0.0
//      for |yDiff| > halfHeight_psf --> return 0.0

// Results from IRAF imexam:
// no shift --> center at 8,8
// xshift=-0.25 --> center at 8.24,8
// xshift=-1 --> center at 9,8
// 
// central pixel is *normally* at 8,8 --> 0,0
// 
// xshift = -1 --> peak-flux center of model now has coordinate 7,8 --> -1,0
// So xshift can be interpreted as "change center of model/PSF from 0,0 to xshift,0"


       
#include <stdio.h>
#include <math.h>
#include <string>
#include "fftw3.h"
#include "gsl/gsl_spline2d.h"

#include "image_io.h"
#include "commandline_parser.h"
#include "utilities_pub.h"
#include "function_objects/function_object.h"
#include "func_pointsource.h"
#include "psf_interpolators.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */

static string  psfFilename1 = "simplegauss.fits";
static string  psfFilename2 = "tests/psf_moffat_35.fits";


typedef struct {
  std::string  inputFilename;
  std::string  outputFilename;
  double  xShift;
  double  yShift;
  double  x0;
  double  y0;
  int  debug;
} commandOptions;





/* ------------------- Function Prototypes ----------------------------- */

// double * MakeShiftedImage2( PsfInterpolator_bicubic *theInterpolator, int nRows, int nColumns, 
// 							double xShift, double yShift, commandOptions *theOptions );
double * MakeShiftedImage3( PointSource *ptSourceFunc, int nRows, int nColumns, 
							double x0, double y0, commandOptions *theOptions );

void ProcessInput( int argc, char *argv[], commandOptions *theOptions );


/* ---------------- MAIN ----------------------------------------------- */

int main( int argc, char *argv[] )
{
  string  psfFilename = psfFilename1;
  string  outputImageFilename = "output_spline2dtest.fits";
  int  nColumns_psf = 0;
  int  nRows_psf = 0;
  long  nPixels_psf = 0;
  double  *psfPixels = nullptr;
  double  *shiftedImage = nullptr;
  commandOptions  options;

  options.inputFilename = psfFilename1;
  options.outputFilename = "output_spline2dtest.fits";
  options.xShift = 0.0;
  options.yShift = 0.0;
  options.x0 = 3.0;
  options.y0 = 3.0;
  options.debug = 0;

  ProcessInput(argc, argv, &options);

//   printf("nCol = 5: ");
//   int nC = 5;
//   double xBound = (nC - 1) / 2.0;
//   for (int n = 0; n < nC; n++)
//     printf("%g, ", n - xBound);
//   printf("\n");
//   printf("nCol = 5: ");
//   nC = 6;
//   xBound = (nC - 1) / 2.0;
//   for (int n = 0; n < nC; n++)
//     printf("%g, ", n - xBound);
//   printf("\n");


  printf("Reading input image (\"%s\") ...\n", options.inputFilename.c_str());
  psfPixels = ReadImageAsVector(options.inputFilename, &nColumns_psf, &nRows_psf);
  if (psfPixels == nullptr) {
    fprintf(stderr,  "\n*** ERROR: Unable to read PSF image file \"%s\"!\n\n", 
    			options.inputFilename.c_str());
    return -1;
  }
  nPixels_psf = (long)(nColumns_psf * nRows_psf);
  printf("naxis1 [# pixels/row] = %d, naxis2 [# pixels/col] = %d; nPixels_tot = %ld\n", 
           nColumns_psf, nRows_psf, nPixels_psf);


  PsfInterpolator *theInterpolator = new PsfInterpolator_bicubic(psfPixels, nColumns_psf, nRows_psf);
  PointSource * pointSourceFunc = new PointSource();
  pointSourceFunc->AddPsfInterpolator(theInterpolator);
  
  // Generate shifted image
  shiftedImage = MakeShiftedImage3(pointSourceFunc, nColumns_psf, nRows_psf, 
  									options.x0, options.y0, &options);

  printf("\nSaving output image (\"%s\") ...\n", options.outputFilename.c_str());
  vector<string>  commentStrings;
  SaveVectorAsImage(shiftedImage, options.outputFilename, nColumns_psf, nRows_psf, commentStrings);

  // Free memory
  delete theInterpolator;
  delete pointSourceFunc;
  if (psfPixels != nullptr)
    fftw_free(psfPixels);
  if (shiftedImage != nullptr)
    fftw_free(shiftedImage);
  
  printf("Done.\n");
}



// double * MakeShiftedImage( double *inputImage, int nColumns, int nRows,
// 							double xShift, double yShift, commandOptions *theOptions )
// {
//   double  *xArray, *yArray;
//   double  *newImage;
//   double  xBound, yBound, deltaXMin, deltaXMax, deltaYMin, deltaYMax;
//   long  nPixelsTot = (long)(nRows * nColumns);
//   int  n, result;
//   size_t xSize = (size_t)(nColumns);
//   size_t ySize = (size_t)(nRows);
// 
//   // construct x and y index arrays (1 .. nColumns or nRows)
//   xBound = (nColumns - 1) / 2.0;
//   yBound = (nRows - 1) / 2.0;
//   xArray = (double *)calloc((size_t)nColumns, sizeof(double));
//   yArray = (double *)calloc((size_t)nRows, sizeof(double));
//   for (n = 0; n < nColumns; n++)
//     xArray[n] = n - xBound;
//   for (n = 0; n < nRows; n++)
//     yArray[n] = n - yBound;
//   deltaXMin = -xBound;
//   deltaXMax = xBound;
//   deltaYMin = -yBound;
//   deltaYMax = yBound;
// 
//   gsl_spline2d * splineInterp;
//   gsl_interp_accel *xacc = gsl_interp_accel_alloc();
//   gsl_interp_accel *yacc = gsl_interp_accel_alloc();
//   splineInterp = gsl_spline2d_alloc(gsl_interp2d_bicubic, xSize, ySize);
//   result = gsl_spline2d_init(splineInterp, xArray, yArray, inputImage, xSize, ySize);
// 
//   newImage = fftw_alloc_real(nPixelsTot);
//   
//   for (long k = 0; k < nPixelsTot; k++) {
//     int  i,j;
//     double  x, y, newVal;
//     j = k % nColumns;
//     i = k / nColumns;
//     y = i - yBound + yShift;
//     x = j - xBound + xShift;
//     if (theOptions->debug > 0)
//       printf("   k = %ld: i,j = %d,%d and x,y = %g,%g\n", k, i, j, x, y);
//     if ((x < deltaXMin) || (x > deltaXMax) || (y < deltaYMin) || (y > deltaYMax))
//       newVal = 0.0;
//     else
//       newVal = gsl_spline2d_eval(splineInterp, x, y, xacc, yacc);
//     newImage[i*nColumns + j] = newVal;
//   }
// 
// 
//   gsl_spline2d_free(splineInterp);
//   gsl_interp_accel_free(xacc);
//   gsl_interp_accel_free(yacc);
//   free(xArray);
//   free(yArray);
//   
//   return newImage;
// }



// double * MakeShiftedImage2( PsfInterpolator_bicubic *theInterpolator, int nRows, int nColumns, 
// 							double xShift, double yShift, commandOptions *theOptions )
// {
//   double  *newImage;
//   long  nPixelsTot = (long)(nRows * nColumns);
//   double  xBound = (nColumns - 1) / 2.0;
//   double  yBound = (nRows - 1) / 2.0;
// 
//   newImage = fftw_alloc_real(nPixelsTot);
// 
// 
//   for (long k = 0; k < nPixelsTot; k++) {
//     int  i,j;
//     double  x, y;
//     j = k % nColumns;
//     i = k / nColumns;
//     y = i - yBound + yShift;
//     x = j - xBound + xShift;
//     if (theOptions->debug > 0)
//       printf("   k = %ld: i,j = %d,%d and x,y = %g,%g\n", k, i, j, x, y);
//     newImage[i*nColumns + j] = theInterpolator->GetValue(x, y);
//   }
//   
//   return newImage;
// }


double * MakeShiftedImage3( PointSource *ptSourceFunc, int nRows, int nColumns, 
							double x0, double y0, commandOptions *theOptions )
{
  double  *newImage;
  long  nPixelsTot = (long)(nRows * nColumns);
  double  paramVect[1] = {1.0};

  newImage = fftw_alloc_real(nPixelsTot);

  ptSourceFunc->Setup(paramVect, 0, x0, y0);

  for (long k = 0; k < nPixelsTot; k++) {
    int  i,j;
    double  x, y;
    j = k % nColumns;
    i = k / nColumns;
    y = i + 1.0;
    x = j + 1.0;
    if (theOptions->debug > 0)
      printf("   k = %ld: i,j = %d,%d and x,y = %g,%g\n", k, i, j, x, y);
    newImage[i*nColumns + j] = ptSourceFunc->GetValue(x, y);
  }
  
  return newImage;
}



void ProcessInput( int argc, char *argv[], commandOptions *theOptions )
{

  CLineParser *optParser = new CLineParser();
  string  tempString = "";

  /* SET THE USAGE/HELP   */
  optParser->AddUsageLine("Usage: ");
  optParser->AddUsageLine("   spline2dtest [options] inputFile [outputFile]");
  optParser->AddUsageLine(" -h  --help                   Prints this help");
//   optParser->AddUsageLine("     --xshift <value>         x-shift in pixels");
//   optParser->AddUsageLine("     --yshift <value>         y-shift in pixels");
  optParser->AddUsageLine("     --x0 <value>             x-center in pixels");
  optParser->AddUsageLine("     --y0 <value>             y-center in pixels");
  optParser->AddUsageLine("");
  optParser->AddUsageLine("     --debug                  extra printouts");

  optParser->AddFlag("help", "h");
  optParser->AddFlag("debug");
//   optParser->AddOption("xshift");
//   optParser->AddOption("yshift");
  optParser->AddOption("x0");
  optParser->AddOption("y0");


  // Comment this out if you want unrecognized (e.g., mis-spelled) flags and options
  // to be ignored only, rather than causing program to exit
  optParser->UnrecognizedAreErrors();
  
  int status = optParser->ParseCommandLine( argc, argv );
  /* Process the results: actual arguments, if any: */
  if (optParser->nArguments() > 0) {
    theOptions->inputFilename = optParser->GetArgument(0);
    printf("\tInput image file = %s\n", theOptions->inputFilename.c_str());
  }
  if (optParser->nArguments() > 1) {
    theOptions->outputFilename = optParser->GetArgument(1);
    printf("\tOutput image file = %s\n", theOptions->outputFilename.c_str());
  }

  /* Process the results: options */
  if (status < 0) {
    printf("\nError on command line... quitting...\n\n");
    delete optParser;
    exit(1);
  }

  /* Process the results: options */
  // First four are options which print useful info and then exit the program
  if ( optParser->FlagSet("help") || optParser->CommandLineEmpty() ) {
    optParser->PrintUsage();
    delete optParser;
    exit(1);
  }
  if ( optParser->FlagSet("debug") || optParser->CommandLineEmpty() ) {
    theOptions->debug = 1;
  }
//   if (optParser->OptionSet("xshift")) {
//     if (NotANumber(optParser->GetTargetString("xshift").c_str(), 0, kAnyReal)) {
//       fprintf(stderr, "*** ERROR: xshift should be a real number!\n");
//       delete optParser;
//       exit(1);
//     }
//     theOptions->xShift = atof(optParser->GetTargetString("xshift").c_str());
//     printf("\txshift = %g pixel\n", theOptions->xShift);
//   }
//   if (optParser->OptionSet("yshift")) {
//     if (NotANumber(optParser->GetTargetString("yshift").c_str(), 0, kAnyReal)) {
//       fprintf(stderr, "*** ERROR: yshift should be a real number!\n");
//       delete optParser;
//       exit(1);
//     }
//     theOptions->yShift = atof(optParser->GetTargetString("yshift").c_str());
//     printf("\tyshift = %g pixel\n", theOptions->yShift);
//   }
  if (optParser->OptionSet("x0")) {
    if (NotANumber(optParser->GetTargetString("x0").c_str(), 0, kAnyReal)) {
      fprintf(stderr, "*** ERROR: x0 should be a real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->x0 = atof(optParser->GetTargetString("x0").c_str());
    printf("\tx0 = %g pixel\n", theOptions->x0);
  }
  if (optParser->OptionSet("y0")) {
    if (NotANumber(optParser->GetTargetString("y0").c_str(), 0, kAnyReal)) {
      fprintf(stderr, "*** ERROR: y0 should be a real number!\n");
      delete optParser;
      exit(1);
    }
    theOptions->y0 = atof(optParser->GetTargetString("y0").c_str());
    printf("\ty0 = %g pixel\n", theOptions->y0);
  }

  delete optParser;

}

