/*! \file
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
void  DownsampleAndReplace( double *oversampledImage, int nOversampCols, int nOversampRows,
							int nOversampPSFCols, int nOversampPSFRows,
							double *mainImage, int nMainCols, int nMainRows,
							int nMainPSFCols, int nMainPSFRows,
							int startX, int startY, int oversampleScale, int debugLevel );

#endif /* _DOWNSAMPLE_H_ */
