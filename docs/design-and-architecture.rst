Design and Architecture of Imfit
================================

General design and operation of imfit
-------------------------------------

The basic operation of imfit, as implemented in imfit\_main.cpp, is: 1.
Process command-line options 2. Read and process configuration file 1.
Optional data-image characteristics (A/D gain, read noise, original-sky
background) 2. Specifications of the model: functions, initial parameter
values and limits 3. Read in user-supplied data from images 1. Image to
be fit 2. Mask image, if any 3. Noise/error/weight image, if any 4. PSF
image(s), if any 4. Create ModelObject instance and supply it with data
1. List of image functions to use 1. The ModelObject instance will then
instantiate corresponding FunctionObject instances 2. PSF image data (if
supplied by user) 3. Image to be fit 4. Data image characteristics 5.
Oversampled PSF image (if supplied by user) 6. Mask image data (if
supplied by user) 7. Type of fit statistic to be calculated 8.
Noise/error/weight image data (if supplied by user) 9. Call FinalSetup()
method on the ModelObject instance 5. Do the fit 1. Set up initial
parameter vector and parameter-limits structure (mp\_par structure) 2.
Call the user-specified solver (Levenberg-Marquardt is default) 7. Print
summary of fit 8. Optionally, do bootstrap resampling to get parameter
uncertainty estimates 9. Save results 1. Save best-fitting parameter
values 2. Optionally, save best-fitting model image and/or residual
image

General design and operation of makeimage
-----------------------------------------

The operation of makeimage is similar to (but simpler than) imfit.

The key difference in terms of generating the model image is that the
user must specify the *size* of the output image, since there is no data
image to use as a reference. There are three ways to do this: 1.
Command-line options (--ncols, --nrows) 2. Specification within the
configuration file (lines beginning NROWS and NCOLS) 3. Use a reference
image (--refimage command-line option)

The basic operation of makeimage, as implemented in makeimage\_main.cpp,
is: 1. Process command-line options 2. Read and process configuration
file 1. Specifications of the model: functions and parameter values 2.
Optional image-size specifications 3. Read in PSF image(s), if supplied
4. Create ModelObject instance and supply it with data 1. Give it list
of image functions to use 2. PSF image data (if supplied by user) 3.
Dimensions of model image 4. Oversampled PSF data (if supplied by user)
5. Create model image: 1. Set up parameter vector 2. Call
CreateModelImage() method of ModelObject instance, passing in the
parameter vector 6. Save the image 1. Generate string vector of comments
for image header 2. Call SaveVectorAsImage() function (image\_io.h).

The ModelObject Class
---------------------

The heart of Imfit is the ModelObject class. This is instantiated once
in imfit (and once in makeimage), and holds (among other things) the
following data: - Model image pixel values - Vector of FunctionObject
instances which define the model - Array of parameter values for the
model - Image pixel data for PSFs, if any - Convolver objects(s) for PSF
convolution - Data image pixel values (when used in imfit) - other
characteristics of data image (size, gain, read noise, etc.) - Error
image (either supplied by user or internally generated) - Mask image, if
supplied by user - Internal weight image (from combination of error
image and weight image).

Not all of these data members are initialized or used -- for example,
none of the data-related members are used by makeimage, and the error
image is not generated or used if a fit uses Poisson-based statistics
instead of chi^2.

The main functionality of the class includes: - Generating a model image
using the vector of FunctionObjects and the current array of parameter
values - Generating individual-function model images (i.e., the
--output-functions option of makeimage) - Computing individual-function
and total fluxes for current model - Computing deviances vector between
current model image and data image (for use by Levenberg-Marquardt
solver) - Computing the fit statistic from comparison of the current
model image and the data image - Printing current parameter values -
Generating a bootstrap resampling pixel-index vector when bootstrap
resampling is being done

Constructing a model image
--------------------------

The actual generation of pixel values in the model image depends on the
vector of FunctionObjects and the current array of corresponding
parameter values. (And, of course, convolution with a PSF if that is
requested.) The FunctionObjects vector contains instantiations of one or
more classes (e.g., Gaussian, Exponential, Sersic, ExponentialDisk3D)
which are subclasses of the abstract base class FunctionObject
(function\_object.h).

When a model image is constructed, the first step is to call the Setup()
method on each FunctionObject and pass in the corresponding parameter
values. This enables the FunctionObject instances to do any initial
computations which don't depend on actual location within the image.

The value of an individual pixel in the model image is obtained by
iterating over the individual FunctionObjects, calling its GetValue()
method with the current pixel coordinates
:raw-latex:`\f`\ :math:`(x,y)\f`, and adding up all the return values.
This all takes place within a loop which iterates over all the pixels in
the model image; this loop is wrapped in OpenMP directives to allow
parallelization across multiple CPU cores.
