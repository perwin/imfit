// Code for setting up ModelObject instances with data

#ifndef _SETUP_MODEL_OBJECT_H_
#define _SETUP_MODEL_OBJECT_H_

#include <vector>
#include <memory>

#include "options_base.h"
#include "model_object.h"
#include "psf_oversampling_info.h"

using namespace std;

static vector<PsfOversamplingInfo *> EMPTY_PSF_OVERSAMPLING_PTR_VECTOR;


ModelObject* SetupModelObject( std::shared_ptr<OptionsBase> options, vector<int> nColumnsRowsVector, 
					double *dataPixels, double *psfPixels=NULL, double *maskPixels=NULL, 
					double *errorPixels=NULL, 
					vector<PsfOversamplingInfo *> psfOversampleInfoVect=EMPTY_PSF_OVERSAMPLING_PTR_VECTOR ); 


#endif  // _SETUP_MODEL_OBJECT_H_
