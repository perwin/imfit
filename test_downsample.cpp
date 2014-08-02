// Code for exercising and testing downsample.cpp
// g++ -o test_downsample test_downsample.cpp downsample.cpp image_io.cpp -I/usr/local/include -lcfitsio -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

using namespace std;

#include "image_io.h"
#include "downsample.h"


// 20x20 image of zeros
string  simpleNullImage_filename = string("simpleimage_20x20_zeros.fits");
// 6x6 image of ones
string  osampOnesImage_filename = string("simpleimage_6x6_ones.fits");
// 6x6 image of ones, with LL corner pixel = 10
string  osampOnesImage3_filename = string("simpleimage_6x6_ones+ten.fits");

// desired output of downsample-and-replace with overSample = 1 -- has [5:11,11:16] = 1.0
string  modifiedNullImage1_filename = string("simpleimage_20x20_zeros+ones6x6.fits");
// desired output of downsample-and-replace with overSample = 3 -- has [5:6,11:13] = 1.0
string  modifiedNullImage2_filename = string("simpleimage_20x20_zeros+ones2x2.fits");
// desired output of downsample-and-replace with overSample = 3, using osampOnesImage3_filename -- has [5:6,11:13] = 1.0,
// except [5,11] = 2.0
string  modifiedNullImage3_filename = string("simpleimage_20x20_zeros+3x3mix.fits");

string  outputImageName = string("bob.fits");



/* ---------------- PUBLIC METHOD: PrintImage ------------------------- */
// Basic function which prints an image to stdout.  Mainly meant to be
// called by PrintInputImage, PrintModelImage, and PrintWeights.

void PrintImage( double *pixelVector, int nColumns, int nRows )
{

  // The following fetches pixels row-by-row, starting with the bottom
  // row (i.e., what we would normally like to think of as the first row)
  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++)   // step by column number = x
      printf(" %f", pixelVector[i*nColumns + j]);
    printf("\n");
  }
  printf("\n");
}



int main(int argc, char *argv[])
{
  int  nColsMain, nRowsMain;
  int  nColsOsamp, nRowsOsamp;
  int  status;
  int  debug = 1;
  int  osampScale = 3;

  double *mainImage = ReadImageAsVector(simpleNullImage_filename, &nColsMain, &nRowsMain);
  double *osampImage = ReadImageAsVector(osampOnesImage3_filename, &nColsOsamp, &nRowsOsamp);

//  printf("The oversampled image:\n");
//  PrintImage(osampImage, nColsOsamp, nRowsOsamp);
  
//  DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,0,0, mainImage, nColsMain,nRowsMain,0,0,
//     					5,11, osampScale, debug);
  DownsampleAndReplace(osampImage, nColsOsamp,nRowsOsamp,0,0, mainImage, nColsMain,nRowsMain,2,2,
     					5,11, osampScale, debug);
  vector<string>  dummy;
  status = SaveVectorAsImage(mainImage, outputImageName, nColsMain, nRowsMain, dummy);

  free(mainImage);
  free(osampImage);
}

