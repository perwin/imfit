The FunctionObject Class and Its Subclasses
===========================================

Overview
--------

Each image function used in the various Imfit programs is a subclass of
an abstract base class called FunctionObject. Subclasses are
defined for specific 2D image functions, each of which computes pixel
intensity values based on input parameters and pixel coordinates.
An instance of the ModelObject class holds a std::vector of pointers
to one or more FunctionObjects, and constructs a model image by
iterating over this vector.

When a model image is being computed, the first step is to call each
FunctionObject instance's Setup() method and pass in the current array
of parameter values, as well as the central position (*x*\ 0,\ *y*\ 0)
for the current function block. In practice, the entire parameter vector
for the model is passed in, along with an integer offset telling the
FunctionObject instance where to look for *its* particular parameter
values.

Then the caller requests intensity values for individual pixels by
repeatedly calling the GetValue() method and passing in the current
pixel coordinates.

The main methods of this class are:

-  Setup() -- Called by ModelObject to pass in the current parameter
   vector at the beginning of the computation of a model image; this
   allows the image function to store the relevant parameter values and
   do any useful computations that don't depend on pixel position.

-  GetValue() -- Called by ModelObject once for each pixel in the model
   image, to pass in the current pixel values (*x*,\ *y*); the image function
   uses these to compute and return the appropriate intensity value for
   that pixel.
   
   Many FunctionObject subclasses have additional private methods that
   do most of the calculations for each GetValue() call.

Note that the base class includes variant virtual methods for possible *1D* 
subclasses, though Imfit does not use any of these.


API
---

**Files:** ``function_objects/function_object.h``, ``function_objects/function_object.cpp``


.. doxygenclass:: FunctionObject
   :members:
   :private-members:
   :protected-members:
