# Getting Started with Imfit


## Fitting Your First Image

Imfit requires, as a bare minimum, two things:
1. An image in FITS format containing the data to be fit;
2. A configuration file describing the model you want to fit to the data.

So to start off, we'll try fitting the image file ic3478rss_256.fits (which is a
256 x 256-pixel cutout from an SDSS r-band image) with a simple exponential
model, which is described in the configuration file config_sersic_ic3478_256.dat.
To do the fit, just type:

imfit ic3478rss_256.fits -c config_exponential_ic3478_256.dat --sky=130.14

(The "--sky-130.14" is a note to imfit that the image has had a background sky
level of 130.14 counts/pixel previously subtracted; if we don't include this,
imfit will get confused by the fact that some of the pixels have slightly negative
values.)

Imfit will print some preliminary information (confirming which files are being
used, the size of the image being fit, the image functions used in the model, etc.).
It will then call the minimization routine, which by default prints a minimal set
of updates for each iteration. At the end, a summary of the fit is printed
(final chi^2, etc.), along with the best-fitting parameters of the model.
These parameters are also saved in a text file: bestfit_parameters_imfit.dat,
while also has a record of how imfit was called and a short summary of the fit.
(You can specify a different name for the output file via the --save-params option.)

Congratulations; you've fit your first image!


Reduced Chi^2 = 0.450366
AIC = 29524.514967, BIC = 29579.055814

X0		128.8530 # +/- 0.0517
Y0		129.1035 # +/- 0.0633
FUNCTION Exponential
PA		19.7364 # +/- 0.467417
ell		0.231428 # +/- 0.00333939
I_0		315.173 # +/- 1.33252
h		20.5726 # +/- 0.0748152


## Inspecting the Fit: Model Images and Residuals

So what kind of 


If you look at the residual image, you can see it's systematically bright in the
center, with an oval region of negative pixels outside. This is a pretty good indication
that the exponential model isn't actually a good match to the data, something we'll
try to address in Section XX.


## Better Fits: Telling Imfit the Truth About the Image

Leaving aside the question of possible mismatch between an exponential model and the
galaxy, this isn't the best possible fit yet for our model.  (You may have noticed
that imfit reported a reduced chi^2 value of ~ 0.45, which a sign
something odd is going on.) For one thing, we've deceived imfit
about the nature of the data. The chi^2 fitting process that imfit uses
is based on the Gaussian approximation to Poisson statistics, and
assumes that the pixel values in the image are detected photoelectrons
(or N-body particles, or something else that obeys Poisson statistics).
In reality, our image deviates from this ideal in three ways:
1. There was a sky background that was subtracted from the image;
2. The pixel values are counts (ADUs), not detected photoelectrons;
3. The image has some Gaussian read noise.

To fix this, we can tell imfit three things:
1. The original background level (which we're already doing, via the --sky option);
2. The A/D gain in electrons/count, via the --gain option;
3. The read noise value (in electrons), via the --readnoise option

In the case of this SDSS image, the corresponding tsField FITS table (from the SDSS
archive) has information about the A/D gain and the read noise (or "dark variance")
and tells us that the gain and read noise are 4.725 and 4.3 electrons, respectively,
for the r-band image.

So we can re-run the fit with the following command:

imfit ic3478rss_256.fits -c config_exponential_ic3478_256.dat --sky=130.14 --gain=4.725 --readnoise=4.3

Now the reduced chi^2 is about 2.1, which isn't necessarily that good, but is at
least statistically plausible!

Reduced Chi^2 = 2.082564
AIC = 136482.400611, BIC = 136536.941458

X0		128.8540 # +/- 0.0239
Y0		129.1028 # +/- 0.0293
FUNCTION Exponential
PA		19.7266 # +/- 0.217212
ell		0.23152 # +/- 0.00155236
I_0		316.313 # +/- 0.619616
h		 20.522 # +/- 0.0346742



## Better Fits: Masking

If you look at the image (e.g., with ds9 or another FITS-displaying
program), you can see features that probably aren't part of the galaxy
-- for example, there are certainly three (and possibly five) distinct,
small objects near the galaxy which are probably foreground stars or
background galaxies. Since they're relatively bright compared to the
outer parts of the galaxy, they will bias the fit.

To prevent this from happening, you can mask out parts of an image. This is
done with a separate mask image: an image of the same size as the data, but with
pixel values = 0 for all the "good" pixels and >= 1 for all the "bad" pixels
(i.e., those pixels you want Imfit to ignore).

The file ic3478rss_256_mask.fits in the examples directory is a mask image. You can
use it in the fit with the "--mask" option:

$ imfit ic3478rss_256.fits -c config_sersic_ic3478_256.dat --mask ic3478rss_256_mask.fits --sky=130.14 --gain=4.725 --readnoise=4.3

The reduced chi^2 is slightly smaller; in addition, the position angle, ellipticity, and

Reduced Chi^2 = 1.964467
AIC = 124602.443320, BIC = 124656.787960

X0		128.8793 # +/- 0.0237
Y0		129.0589 # +/- 0.0289
FUNCTION Exponential
PA		18.7492 # +/- 0.23086
ell		0.220646 # +/- 0.00159077
I_0		321.631 # +/- 0.634224
h		20.0684 # +/- 0.034584


## Better Fits: Trying a different model



 

## Better Fits: PSF Convolution



## Chi-Squared and All That: Using Different Fit Statistics

Chi-squared from data, chi-squared from model

Poisson MLR


# Bits of Advice
