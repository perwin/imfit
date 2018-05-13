Convolver
=========

Overview
--------

This class holds data and functions relating to FFT-based convolution of
a point-spread-function (PSF) image with an input image (e.g., the model
image computed by ModelObject).

An instance of this class is set up and used by main's ModelObject instance
to handle standard PSF convolution of the whole model image. In addition,
if there are any oversampled regions, then each OversampledRegion instance
will have its own Convolver instance to handle convolution with the
appropriate (oversampled) PSF image.


Setup and Use
-------------

To set up a Convolver object, you first call its SetupPSF method and
pass in the PSF image data, size, and whether the PSF should be
normalized. Then you call the SetupImage method to tell the Convolver
object about the size of the image that will be convolved with the PSF,
and finally you call the DoFullSetup to tell the Convolver object to
do the necessary allocations and FFT setup.


To use the Convolver object, you simply call its ConvolveImage method
with the image array as input; the input image will be updated in place.

(In Imfit, this is all done within the ModelObject class.)


API
---

**Files:** ``core/convolver.h``, ``core/convolver.cpp``


.. doxygenclass:: Convolver
   :members:
   :private-members:
   :protected-members:
