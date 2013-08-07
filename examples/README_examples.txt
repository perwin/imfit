FILES:
The files in this directory are sample images and configuration files to give
you a quick example of how to use Imfit.

* Galaxy image (a 256x256-pixel cutout from an SDSS DR7 r-band image of the dE galaxy
IC 3478): ic3478rss_256.fits
(this image has had a constant sky value of 130.14 subtracted from it)

* Mask file: ic3478rss_256_mask.fits
(the DS9 region file which was used to generate the mask is ic3478rss_256.reg)

* PSF image: psf_moffat_51.fits
(a 51x51-pixel image)

* Imfit configuration file for fitting a single elliptical Sersic function:
config_sersic_ic3478_256.dat


IMFIT QUICKSTART:
To fit the galaxy image using a single Sersic function, without any PSF convolution,
but using the mask (and SDSS gain, read noise, and original sky background value
to estimate the noise image):

$ imfit ic3478rss_256.fits -c config_sersic_ic3478_256.dat --mask ic3478rss_256_mask.fits \
--gain=4.725 --readnoise=4.3 --sky=130.14

This converges to a fit in about 2.5 seconds on a 2009 MacBook Pro (2.8
GHz Core 2 Duo processor), or less than a second on a more modern
machine. The fit will be saved in a text file called (by default)
bestfit_parameters_imfit.dat


Same, but now using PSF convolution:

$ imfit ic3478rss_256.fits -c config_sersic_ic3478_256.dat --mask ic3478rss_256_mask.fits \
--gain=4.725 --readnoise=4.3 --sky=130.14 --psf psf_moffat_51.fits


The PSF image is a simple Moffat function (matched to stars in the SDSS
image); it was generating using the makeimage program and the
configuration file config_makeimage_moffat_psf.dat (which specifies an
output image 51 x 51 pixels on a side)

$ makeimage config_makeimage_moffat_psf.dat -o psf_moffat_51.fits

