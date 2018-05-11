General Notes and Advice for Fitting Images
===========================================

**WARNING: this page is currently an incomplete work in progress!**

Pixel values and uncertainties
------------------------------

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
tell imfit about the error image (via the ``--noise`` command-line
option, plus the ``--errors-are-variances`` flag if the pixel values in
the error image are variances instead of sigmas).

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order for imfit to correctly estimate the per-pixel uncertainties, it
is *very important* that the pixel values be convertible to units of the
original detected quantities -- i.e., integrated photo-electrons for a
typical detector, or particles/pixel for an *N*-body simulation.

Since the actual recording of an astronomical image almost always
involves conversion of detected photo-electrons to counts (also called
"data numbers" or ADUs) via an A/D converter, the output image is
already no longer in

So you must tell imfit how to convert the pixel values back to the
original values, at least in a general, image-wide sense. (It's probably
overkill to worry about the effects of flat-fielding.)

The simple way to do this is to give imfit an effective gain factor, via
the GAIN parameter in the config file, or the ``--gain`` command-line
option. The effective gain is whatever number will convert pixel values
back to photo-electrons. For a processed (but not photometrically
calibrated) image, where the pixel values are in ADUs, this is simply
the gain of the original A/D conversion, in units of electrons/ADU.
Imfit will then multiply the pixel values by the GAIN parameter to turn
them back into photo-electrons. (Note that some data-taking setups will
record an inverse "gain" in units of ADU/electrons; you will need to
invert this before passing it to imfit!)

If your image has photometrically calibrated flux units (Jy,
nanomaggies, whatever), then the GAIN parameter must be a number that
will convert these back to photo-electrons. XXX If an image happens to
be in units of magnitudes per square arc second, then it must be
converted to *linear* intensity values.

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
   ``ORIGINAL_SKY`` (or ``--sky`` commandline option) parameter for
   imfit;

-  OR: Subtract 1000 from the image, then include a background component
   (e.g., FlatSky) in the model with initial value = 210.7 (possibly
   fixed to that value, if you don't want imfit to vary the background).

The A/D gain values are included in the tsField FITS tables that go with
each field; typical values can also be found in the table at the bottom
of `this
page <http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html>`__

DR 8+ Images
~~~~~~~~~~~~

Images that are part of DR8 and later releases are more problematic,
because they have been preprocessed in slightly complicated and
potentially confusing ways.

First, the pixel values have been photometrically calibrated and
transformed into flux values -- specifically, "nanomaggies" (units of
3.631e-6 Jy). To fit the images properly, you will need to convert these
back to counts, or else combine the flux conversion with the A/D gain to
make an "effective gain" parameter as input to imfit.

The photometric calibration is recorded in the image header using the
NMGY header keyword. You can convert the image back to ADUs via

::

    image_ss_counts = image_ss_nmgy / NMGY

Second, a 2D model of the background is computed and subtracted from the
image as part of the pipeline, so you will definitely need to include
some approximation of the background value as an ``ORIGINAL_SKY``
parameter in the configuration file (or via the ``--sky`` command-line
option). The 2D background model is included in the image files in a
somewhat obscure fashion, as a series of floating-point values for
interpolation, stored as a table in extension 2 of the An additional
problem is that these sky-model values are in counts, *not* in the
nanomaggy units of the processed data image.

An estimate of the original sky background can be obtained by averaging
the interpolation values in the second extension. For example, using
``numpy`` and the ``astropy.io.fits`` module in Python:

::

    >>> import numpy as np
    >>> from astropy.io import fits
    >>> hdu_list = fits.open('<path-to-SDSS-image-file>')
    >>> sky_bintable = hdu_list[2]
    >>> np.mean(sky_bintable.data['ALLSKY'])

The simplest approach is probably to convert the image from nanomaggies
to counts, then use the standard A/D gain values and the mean sky value
from the second FITS extension.

Alternately, you could combine the NMGY value with the gain to get an
effective-gain value that converts the nanomaggie values directly to
photo-electrons ("gain" = A/D gain / NMGY) -- but then you will have to
convert the sky value from counts to nanomaggies.
