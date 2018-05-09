The ModelObject Class
=====================

Overview
--------

The heart of Imfit is the ModelObject class. This holds information
about an image model and can compute images using the model; it can also
hold a data image and its associated information, and can compute chi^2
and and other statistics by comparing the model and data images.

This class is instantiated once in each program, and holds (among other
things) the following objects and data ("**imfit only**\ " denotes data
used in imfit or imfit-mcmc, but not needed when used in makeimage)

-  Model image pixel values

-  Vector of FunctionObject instances which define the model

-  Array of parameter values for the model

-  Image pixel data for PSFs, if any

-  Convolver objects(s) for PSF convolution

-  **imfit only**: Data image pixel values

-  **imfit only**: other characteristics of data image (size, gain, read
   noise, etc.)

-  **imfit only**: Error image (either supplied by user or internally
   generated)

-  **imfit only**: Mask image, if supplied by user

-  **imfit only**: Internal weight image (from combination of error
   image and weight image).

Not all of these data members are initialized or used -- for example,
none of the data-related ("**imfit only**\ ") members are used by
makeimage, and the error image is not generated or used if a fit uses
Poisson-based statistics instead of chi^2.

The main functionality of the class includes:

-  Generating a model image using the vector of FunctionObjects and the
   current array of parameter values

-  Generating individual-function model images (i.e., the
   --output-functions option of makeimage)

-  Computing individual-function and total fluxes for current model

-  Computing deviances vector between current model image and data image
   (for use by Levenberg-Marquardt solver)

-  Computing the fit statistic (chi^2, etc.) from comparison of the
   current model image and the data image

-  Printing current parameter values

-  Generating a bootstrap resampling pixel-index vector when bootstrap
   resampling is being done



API
---

[**Warning**: This is currently somewhat incomplete!]

**Files:** ``core/model_object.h``, ``core/model_object.cpp``


.. doxygenclass:: ModelObject
   :members:
   :private-members:
   :protected-members:
