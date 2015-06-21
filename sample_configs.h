/** @file
    \brief Header file containing text for example config files, and functions
    for writing them to files.
 *
 */

#ifndef _SAMPLE_CONFIGS_H_
#define _SAMPLE_CONFIGS_H_

#include <cstdio>
#include <string>

std::string configImfitFile = "config_imfit_sample.dat";
std::string configMakeimageFile = "config_makeimage_sample.dat";

/// Full text for the example imfit config file
// compiler will concatenate all this together into a single string
std::string configImfitText = "# This is a sample configuration file for imfit (it will also work with makeimage).\n"
"# Comments are introduced with the \"#\" character; anything following this character\n"
"# on the same line is ignored by imfit. Blank lines are also ignored.\n"
"# Note that to use this file with makeimage, the desired dimensions of the output\n"
"# model image will also need to be supplied.\n\n"
"# Optional image-description parameters for imfit and their default values;\n"
"# uncomment and edit if needed\n"
"# GAIN            1.0   # A/D gain in e-/ADU\n"
"# READNOISE       0.0   # read noise in electrons\n"
"# ORIGINAL_SKY    0.0   # original background value (ADU/pixel) previously subtracted\n"
"# EXPTIME         1.0   # only needed if image has ADU/sec/pixel instead of ADU/pixel\n"
"# NCOMBINED       1     # only needed if image is average of multiple images\n\n\n"
"# A single function block, containing a single image function;\n"
"# examples of optional lower and upper limits for each parameter value are in the 3rd column\n"
"X0   129.0    125,135     # pixel x-coordinate of function center\n"
"Y0   129.0    125,135     # pixel y-coordinate of function center\n"
"FUNCTION Sersic\n"
"PA    18.0    0,90        # position angle in degrees (CCW from +y axis of image)\n"
"ell    0.2    0,1         # ellipticity of isophotes\n"
"n      1.5    0,5         # Sersic index\n"
"I_e    15     0,500       # intensity (ADU/pixel) at r_e\n"
"r_e    25     0,100       # half-light radius (pixels)\n";

/// Full text for the example makeimage config file
std::string configMakeimageText = "# This is a sample configuration file for makeimage.\n"
"# Comments are introduced with the \"#\" character; anything following this character\n"
"# on the same line is ignored by makeimage. Blank lines are also ignored.\n"
"# (Note that this file can also be used with imfit.)\n\n"
"# Image-description parameters for makeimage (image dimensions)\n"
"NCOLS   200      # number of columns (width) in output image\n"
"NROWS   200      # number of rows (height) in output image\n\n\n"
"# A single function block, containing a single image function\n"
"X0   129.0       # pixel x-coordinate of function center\n"
"Y0   129.0       # pixel y-coordinate of function center\n"
"FUNCTION Sersic\n"
"PA    18.0       # position angle in degrees (CCW from +y axis of image)\n"
"ell    0.2       # ellipticity of isophotes\n"
"n      1.5       # Sersic index\n"
"I_e    15        # intensity (ADU/pixel) at r_e\n"
"r_e    25        # half-light radius (pixels)\n";


/// Saves an example config file for imfit. Returns 0 for success, -1 if file couldn't be opened.
int SaveExampleImfitConfig( )
{
  FILE *outFile_ptr = fopen(configImfitFile.c_str(), "w");
  if (outFile_ptr) {
    fprintf(outFile_ptr, "%s", configImfitText.c_str());
    fclose(outFile_ptr);
    return 0;
  }
  else {
    fprintf(stderr, "*** Unable to save sample config file! ***\n");
    return -1;
  }
}


/// Saves an example config file for makeimage. Returns 0 for success, -1 if file couldn't be opened.
int SaveExampleMakeimageConfig( )
{
  FILE *outFile_ptr = fopen(configMakeimageFile.c_str(), "w");
  if (outFile_ptr) {
    fprintf(outFile_ptr, "%s", configMakeimageText.c_str());
    fclose(outFile_ptr);
    return 0;
  }
  else {
    fprintf(stderr, "*** Unable to save sample config file! ***\n");
    return -1;
  }
}


#endif    // _SAMPLE_CONFIGS_H_
