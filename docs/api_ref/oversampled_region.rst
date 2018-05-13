OversampledRegion
=================

Overview
--------

This class handles computation of oversampled regions in the model
image; there will be one instance for each oversampled region requested
by the user (these are stored in ModelObject as a vector of pointers to
individual OversampledRegion objects).

It holds information about the location, size, and pixel
oversampling factor of the oversampled region; it also holds a Convolver
object which performs PSF convolution using the appropriate oversampled
PSF image.

Code contained in this class can generate an oversampled model-image
region, using the FunctionObjects vector passed in by the caller (i.e.,
main's ModelObject instance); this subimage is then convolved with the
oversampled PSF and finally block-downsampled back to the main image's
pixel scale and copied into the corresponding pixels of the main model
image. OpenMP compiler directives are used to speed up computation of
the oversampled model subimage.


Setup and Use
-------------

To set up an OversampledRegion object, you first call its AddPSFVector
method to pass in the PSF image data, size, and whether the PSF should
be normalized. Then you tell the OversampledRegion object about the size
of the main image, the location of the oversampled region, and the
oversampling scale by calling the SetupModelImage method.

(See ModelObject::AddOversampledPSFVector for an example of how this is done
in practice.)

To *use* an OversampledRegion object in the model-image-generation
process, you call its ComputeRegionAndDownsample method, passing in a
pointer to the current model image and the vector of FunctionObjects
which describe the model. You should be sure to update the
FunctionObjects with the current set of parameter values before doing
this. (See ModelObject::CreateModelImage for an example of use.)

(In Imfit, this is all done within the ModelObject class.)



API
---

**Files:** ``core/oversampled_region.h``, ``core/oversampled_region.cpp``


.. doxygenclass:: OversampledRegion
   :members:
   :private-members:
   :protected-members:
