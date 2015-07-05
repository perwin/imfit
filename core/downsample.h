/** @file
    \brief Utility function taking an oversampled sub-region, downsampling it to
           parent image's pixel scale, and copying it into parent image.
 *
 */
/*    Utility functions taking an oversampled sub-region, downsampling it to
 * parent image's pixel scale, and copying it into parent image.
 *
 */

#ifndef _DOWNSAMPLE_H_
#define _DOWNSAMPLE_H_


/// \brief Replaces subsection of main region with oversampled version of subsection,
///        downsampled to match main image scale
void DownsampleAndReplace( double *oversampledImage, const int nOversampCols, 
						const int nOversampRows, const int nOversampPSFCols, 
						const int nOversampPSFRows,	double *mainImage, const int nMainCols, 
						const int nMainRows, const int nMainPSFCols, const int nMainPSFRows,
						const int startX, const int startY, const int oversampleScale, 
						const int debugLevel );

#endif /* _DOWNSAMPLE_H_ */
