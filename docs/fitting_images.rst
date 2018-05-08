General Notes and Advice for Fitting Images
===========================================

It is important to have some idea of the underlying statistical model
for your data.

The basic assumption is that the pixel values in a data image are the
outcome of some statistical sampling process (usually Poisson or
Gaussian) operating on an underlying model of the intensity
distribution.

A very simple example of an underlying model might be a constant sky
background plus an elliptical 2D Sersic function representing a galaxy.
Sampling of this by a 2D detector results in pixel values which are a
combination of the image model and a sampling model for the distribution
of individual intensity value recorded for each pixel.

If your data processing pipeline has produced some kind of error image
(sigma or variance values) which you trust, then you can go ahead and
tell imfit about the error image (via the --noise command-line option,
plus the --errors-are-variances flag if the pixel values in the error
image are variances instead of sigmas).

For the common case where you do *not* have an error image, you can have
imfit estimate the per-pixel uncertainties for you, either from the data
image or from the model. The question is then whether to use proper
Poisson statistics or the usual Gaussian approximation of Poisson
statistics.

If the original count levels (including background) were high (say,
greater than 100 photo-electrons per pixel), then you can probably
assume the Gaussian approximation of Poisson statistics is OK, and use
some variant of chi^2 as the fit statistic.

IMPORTANT NOTE: Converting your pixel values to photo-electrons
---------------------------------------------------------------

In order for imfit to correctly estimate the per-pixel uncertainties, it
is *very important* that the pixel values be in units of the original
detected quantities -- i.e., integrated photo-electrons for a typical
detector, or particles/pixel for an *N*-body simulation.

So you must tell imfit how to convert the pixel values in the data image
from their current values back to the original values, at least in a
general, image-wide sense. (It's probably overkill to worry about the
effects of flat-fielding.)

The simple way to do this is to give imfit an effective gain factor, via
the GAIN parameter in the config file, or the --gain command-line
option. The effective gain is whatever number will convert pixel values
back to photo-electrons. For a processed (but not photometrically
calibrated) image, where the pixel values are in ADUs, gain has units of
electrons/ADU. In typical image processing, the gain has units of
electrons/ADU (so the recorded image has units of ) XXX

If your image has photometrically calibrated flux units (Jy,
nanomaggies, whatever), then the gain parameter must be a number that
will convert these back to photo-electrons. XXX

Special Notes for SDSS Images
-----------------------------

DR7 (and earlier) Images
~~~~~~~~~~~~~~~~~~~~~~~~

DR7 images are fairly straightforward; the only real potential for
confusion is with the background level.

DR7 images do *not* have background subtraction applied (a crude
estimate of the background level is included in the header of each
image). However, each image has an artifical "soft bias" level of 1000
added to each pixel. Thus, if an image has a mean background level of
1210 counts/pixel, the *actual* observed sky background is only 210
counts/pixel. The additional constant value of 1000 should be removed
*and not included in any sky-background levels input to imfit*.

Thus, if you determine that a particular image has a mean background
level of 1210.7 counts/pixel, you can either:

-  Subtract 1210.7 from the image, and use 201.7 (*not* 1210.7) as the
   ``ORIGINAL_SKY`` (or ``--original-sky`` commandline option) parameter
   for imfit

-  Subtract 1000 from the image, then include a background component
   (e.g., FlatSky) in the model with initial value = 210.7 (possibly
   fixed to that value, if you don't want imfit to vary the background).

The A/D gain values are included in the tsField FITS tables that go with
each field.

DR 8+ Images
------------

Images that are part of DR8 and later releases are more problematic,
because they have been preprocessed in slightly complicated and
potentially confusing ways.

First, the pixel values have been photometrically calibrated and
transformed into flux values -- specifically, "nanomaggies" (units of
XXX Jy). To fit the images properly, you will need to convert these back
to counts/pixel, or else combine the flux conversion with the A/D gain
to make an "effective gain" parameter as input to imfit.

The photometric calibration is recorded in the image header using the
NMGY header keyword. You can convert the image back to ADUs xxx

::

    image_ss_adu = image_ss_nmgy / NMGY

Second, a 2D model of the background is computed and subtracted from the
image as part of the pipeline. The 2D background model is included in
the image files in a somewhat obscure fashion: XXX. An additional
problem is that these sky-model values are in ADUs (counts/pixel), *not*
in the nanomaggy units of the processed data image.

The choices are basically:

1. Convert the image from nanomaggies to counts/pixel and use the
   standard A/D gain values + use the mean sky value from the XXX FITS
   extension

2. Compute an effective gain as XXX + convert the mean sky value from
   counts/pixel to nanomaggies
