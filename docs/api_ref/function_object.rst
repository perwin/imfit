The FunctionObject Class and Its Subclasses
===========================================

Overview
--------

Each image function used in the various Imfit programs is a subclass of
the abstract base class FunctionObject. Subclasses of this class are
defined for specific 2D image functions, each of which computes pixel
intensity values based on input parameters and pixel coordinates.
An instance of the ModelObject class holds a std::vector of pointers
to one or more FunctionObjects, and constructs a model image by
iterating over this vector.

The general use when a model image is being computed is to first call a
FunctionObject instance's Setup() method and pass in the current array
of parameter values (including the central position (*x0*,*y0*) for the
current function block). Then the caller (e.g., one of the methods in
ModelObject) requests intensity values for individual pixels by
repeatedly calling the GetValue() method and passing in the current
pixel coordinates.

(In practice, ModelObject first iterates over its internal vector of
FunctionObjects, calling the Setup() method on each; then steps through
the pixels of the model image, calling the GetValue() methods for each
FunctionObject and adding up their values to get the final pixel values.)

The main methods of this class are:

-  Setup() -- Called by ModelObject to pass in the current parameter
   vector at the beginning of the computation of a model image; this
   allows the image function to store the relevant parameter values and
   do any useful computations that don't depend on pixel position.

-  GetValue() -- Called by ModelObject once for each pixel in the model
   image, to pass in the current pixel values (x,y); the image function
   uses these to compute and return the appropriate intensity value for
   that pixel.
   
   Many FunctionObject subclasses have additional private methods that
   do most of the calculations for each GetValue() call.

Note that the base class include variant virtual methods for possible *1D* 
subclasses, though Imfit does not use any of these.


API
---

**Files:** ``function_objects/function_object.h``, ``function_objects/function_object.cpp``


.. doxygenclass:: FunctionObject
   :members:
   :private-members:
   :protected-members:
