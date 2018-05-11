The PsfOversamplingInfo Class
=============================

This is a data-container class meant to assist in the setup phase of
imfit, makeimage, etc. It holds information about a single
user-specified oversampling region: the location and of the region
within the image, its size, the desired pixel oversampling scale, and
the oversampled PSF image to be used with that region.

In practice, a vector of one or more PsfOversamplingInfo objects is set
up in main() based on the user inputs and is then passed to the
SetupModelObject function, where the individual PsfOversamplingInfo
objects are passed to the ModelObject instance via its
AddOversampledPsfInfo method.


API
---

**Files:** ``core/psf_oversampling_info.h``, ``core/psf_oversampling_info.cpp``


.. doxygenclass:: PsfOversamplingInfo
   :members:
   :private-members:
   :protected-members:
