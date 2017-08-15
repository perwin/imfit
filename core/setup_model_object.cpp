// Code for setting up ModelObject instances with data

// Copyright 2017 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <stdlib.h>

#include "setup_model_object.h"
#include "options_base.h"
#include "model_object.h"

using namespace std;

// We assume that the input nColumnsRowsVector has the following entries:
// nColumnsRowsVector[0] = nColumns
// nColumnsRowsVector[1] = nRows
// nColumnsRowsVector[2] = nColumns_psf  [optional]
// nColumnsRowsVector[3] = nRows_psf  [optional]

ModelObject* SetupModelObject( OptionsBase *options, vector<int> nColumnsRowsVector, 
					double *dataPixels, double *psfPixels, double *maskPixels, 
					double *errorPixels, vector<PsfOversamplingInfo *> psfOversampleInfoVect )
{
  ModelObject *newModelObj;
  int  status;
  bool  doingFit_or_MCMC = false;
  int  nColumns, nRows, nColumns_psf, nRows_psf;
  long  nPixels_data, nPixels_psf;
  
  newModelObj = new ModelObject();
  if (options->debugLevel > 0)
    printf("\nSetupModelObject: Setting ModelObject debug level to %d\n\n", options->debugLevel);
  newModelObj->SetDebugLevel(options->debugLevel);

  if (options->maxThreadsSet)
    newModelObj->SetMaxThreads(options->maxThreads);
  newModelObj->SetDebugLevel(options->debugLevel);


  // Add PSF image vector, if present (needs to be added prior to image data or
  // model-image setup, so that ModelObject can figure out proper internal model-image 
  // size when we call SetupModelImage or AddImageDataVector)
  if (options->psfImagePresent) {
    nColumns_psf = nColumnsRowsVector[2];
    nRows_psf = nColumnsRowsVector[3];
    nPixels_psf = nColumns_psf * nRows_psf;
    status = newModelObj->AddPSFVector(nPixels_psf, nColumns_psf, nRows_psf, psfPixels,
    									options->normalizePSF);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddPSFVector!\n\n");
  	  exit(-1);
    }
  }

  nColumns = nColumnsRowsVector[0];
  nRows = nColumnsRowsVector[1];
  nPixels_data = (long)nColumns * (long)nRows;
  if (dataPixels == NULL) {
    // No data image, so we're in makeimage mode
    status = newModelObj->SetupModelImage(nColumns, nRows);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::SetupModelImage!\n\n");
      exit(-1);
    }
  } else {
    // data image exists, so we're in imfit or imfit-mcmc mode
    doingFit_or_MCMC = true;
    status = newModelObj->AddImageDataVector(dataPixels, nColumns, nRows);
    if (status < 0) {
      // Possible failure if attempt to allocate memory for model image fails
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddImageDataVector!\n\n");
      exit(-1);
    }
    newModelObj->AddImageCharacteristics(options->gain, options->readNoise, options->expTime, 
    						options->nCombined, options->originalSky);
  }

  // Add oversampled PSF image vector(s) and corresponding info, if present
  if (options->psfOversampling) {
    for (int i = 0; i < (int)psfOversampleInfoVect.size(); i++) {
      status = newModelObj->AddOversampledPsfInfo(psfOversampleInfoVect[i]);
      if (status < 0) {
        fprintf(stderr, "*** ERROR: Failure in ModelObject::AddOversampledPSFVector!\n\n");
  	    exit(-1);
      }
    }
  }

  // If user supplied a mask image, add it and apply it to the internal weight image
  if (options->maskImagePresent) {
    status = newModelObj->AddMaskVector(nPixels_data, nColumns, nRows, maskPixels,
                             options->maskFormat);
    if (status < 0) {
      fprintf(stderr, "*** ERROR: Failure in ModelObject::AddMaskVector!\n\n");
  	  exit(-1);
    }
  }
  
  
  if (doingFit_or_MCMC) {
    // Specify which fit statistic we'll use, and add user-supplied error image if
    // it exists and we're using chi^2; also catch special case of standard Cash
    // statistic + L-M minimizer
    if (options->useCashStatistic) {
      if ((options->solver == MPFIT_SOLVER) && (! options->printFitStatisticOnly)) {
        fprintf(stderr, "*** ERROR -- Cash statistic cannot be used with L-M solver!\n\n");
        exit(-1);
      }
      status = newModelObj->UseCashStatistic();
      if (status < 0) {
        fprintf(stderr, "*** ERROR: Failure in ModelObject::UseCashStatistic!\n\n");
        exit(-1);
      }
    } 
    else if (options->usePoissonMLR) {
      newModelObj->UsePoissonMLR();
    }
    else {
      // normal chi^2 statistics, so we either add error/noise image, or calculate it
      if (options->noiseImagePresent)
        newModelObj->AddErrorVector(nPixels_data, nColumns, nRows, errorPixels,
                               options->errorType);
      else {
        if (options->useModelForErrors) {
          printf("* No noise image supplied ... will generate noise image from model image.\n");
          status = newModelObj->UseModelErrors();
          if (status < 0) {
            fprintf(stderr, "*** ERROR: Failure in ModelObject::UseModelErrors!\n\n");
            exit(-1);
          }
        }
        else {
          // default mode
          printf("* No noise image supplied ... will generate noise image from input data image.\n");
        }
      }
    }
  }


  return newModelObj;
}
