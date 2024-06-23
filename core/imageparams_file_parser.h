/** @file
    \brief Public interfaces for code which parses and re-arranges image-description
    parameters (pixelScale, rotation, intensityScale, etc.) from image-info files

 */

#ifndef _IMAGEPARAMS_FILE_PARSER_H_
#define _IMAGEPARAMS_FILE_PARSER_H_

#include <string>
#include <vector>
#include <tuple>

#include "definitions_multimage.h"
#include "param_struct.h"

using namespace std;

// Error codes returned by VetInputLines
const int IMAGEPARAMS_FILE_ERROR_NOLINES  = -1;
const int IMAGEPARAMS_FILE_ERROR_NOLINENUMBERS  = -2;
const int IMAGEPARAMS_FILE_ERROR_NOIMAGESTART = -3;
const int  IMAGEPARAMS_FILE_ERROR_BADLINE = -4;


bool VetInputLines( const vector<string> &lines, const vector<int> &lineNumbers,
					vector<int> &errorCodes, vector<int> &badLineNumbers );

std::tuple<double, mp_par, int> GetParameterAndLimits( const string theLine );

// two versions of this: first is simpler version for use by makemultimates
// (no parameter limit info), second is for use by multimfit (includes parsing
// and returning of parameter limit info)
// Returns tuple of (array of image-description parameter values for all images, status)
std::tuple<double *, int>  GetImageDescriptionParams( const vector<ImageInfo> imageInfoList );

std::tuple<double *, int>  GetImageDescriptionParams( const vector<ImageInfo> imageInfoList,
											vector<mp_par>& imageParamLimits,
											bool& paramLimitsExist );

int StoreInfoFromLine( const string theLine, ImageInfo& imageInfo );

int AddInfoFromLines( const vector<string>& lines, vector<ImageInfo>& imageInfoList );

int ReadImageParamsFile( const string& configFileName, int& nImages,
						vector<ImageInfo>& imageInfoList );


#endif  // _IMAGEPARAMS_FILE_PARSER_H_
