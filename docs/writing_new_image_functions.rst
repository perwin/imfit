How to Write New Image Functions
================================

Overview
--------

Imfit comes with a variety of 2D surface-brightness functions ("image
functions"). But it is designed to make it relatively simple to add
*new* image functions, to describe objects or substructures which are
not well modeled by the existing functions. This is done by writing
additional C++ code and re-compiling the various Imfit programs
(``imfit``, ``imfit-mcmc``, ``makeimage``) to include the new function
or functions.

This document provides some guidelines on how to write new image
functions (see also the "Rolling Your Own Image Functions" section of
the Imfit manual, from which this document is excerpted). Even if you're
not very familiar with C++, it should hopefully not be too difficult to
write a new image function, since it can be done by copying and
modifying one of the existing image functions.

A new image function is written in C++ as a subclass of the
FunctionObject base class, which is declared and defined in the
``function_object.h`` and ``function_object.cpp`` source-code files (in
the ``function_objects/`` subdirectory of the source-code distribution).

A new image-function class should provide its own implementation of
(i.e., "override") the following public FunctionObject methods, which
are defined as virtual methods in the base class:

-  The class constructor -- in most cases, the code for this can be
   copied from any of the existing FunctionObject subclasses, unless
   some special extra initialization is needed.
-  ``Setup()`` -- this is used by the calling program to supply the
   current set of function parameters (including the (*x*\ 0,\ *y*\ 0)
   pixel values for the center) prior to determining intensity values
   for individual pixels. This is a convenient place to do any general
   calculations which don't depend on the exact pixel (*x*,\ *y*)
   values.
-  ``GetValue()`` -- this is used by the calling program to obtain the
   surface brightness for a given pixel location (*x*,\ *y*). In
   existing FunctionObject subclasses, this method often calls other
   (private) methods to handle details of the calculation.
-  ``GetClassShortName()`` -- this is a class function which returns the
   short version of the class name as a string.

The new class should also redefine the following internal class
constants:

-  ``N_PARAMS`` --- the number of input parameters (*excluding* the
   central pixel coordinates);
-  ``PARAM_LABELS`` --- a vector of string labels for the input
   parameters;
-  ``FUNCTION_NAME`` --- a short string describing the function;
-  ``className`` --- a string (no spaces allowed) giving the official
   name of the function.

The ``add_functions.cpp`` file should then be updated by:

-  including the header file for the new class;
-  adding 2 lines to the ``PopulateFactoryMap()`` function to add the
   ability to create an instance of the new class.

Finally, the name of the C++ implementation file for the new class
should be added to the ``SConstruct`` file to ensure it gets included in
the compilation; the easiest thing is to add the file's name (without
the ``.cpp`` suffix) to the multi-line ``functionobject_obj_string``
string definition.

Existing examples of FunctionObject subclasses can be found in the
``function_objects/`` subdirectory of the source-code distribution, and
are the best place to look in order to get a better sense of how to
implement new FunctionObject subclasses.

A Simple Example
----------------

To demonstrate the basics of writing a new image function, we'll modify
the existing Gaussian class to make a class called NewMoffat, which
produces an elliptical structure with a Moffat radial profile. (Such a
function already exists in Imfit under the name Moffat, so this is
really a redundant exercise.)

Three basic changes to the existing ``func_gauss.h/cpp`` files are
needed:

1. Change the class name (in this case, from "Gaussian" to "NewMoffat");

2. Change the code which actually computes the function;

3. Add, rename, and delete class data members to accommodate the new
   algorithm.

1. Create and Edit the Header File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming we're in the Imfit source code directory, ``cd`` to the
``function_objects/`` subdirectory, copy the header file
``func_gaussian.h`` to ``func_new-moffat.h``, and end the file to change
the following lines:

::

    #define CLASS_SHORT_NAME "Gaussian"

    class Gaussian : public FunctionObject

    Gaussian( );

to

::

    #define CLASS_SHORT_NAME "NewMoffat"

    class NewMoffat : public FunctionObject

    NewMoffat( );

2. Create and Edit the Class File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy the file ``func_gaussian.cpp`` to ``func_new-moffat.cpp``.

A. Change the following lines in the beginning of the file:

::

    #include "func_gaussian.h"

    const in N_PARAMS = 4;
    const char PARAM_LABELS[][20] = {"PA", "ell", "I_0", "sigma"};
    const char FUNCTION_NAME[] = "Gaussian function";

to

::

    #include "func_new-moffat.h"

    const in N_PARAMS = 5;
    const char PARAM_LABELS[][20] = {"PA", "ell", "I_0", "fwhm", "beta"};
    const char FUNCTION_NAME[] = "Moffat function";

B. In the remainder of the file, change all references to the class name
from ``Gaussian`` to ``NewMoffat`` (e.g., ``Gaussian::Setup`` becomes
``NewMoffat::Setup``).

C. Change the ``Setup`` method. Here, you'll need to change how the
input parameter vector is converted into individual parameters, and do
any useful pre-computations (i.e., computations that depend on the
parameter values, but not on individual pixel values or values derived
from the pixel values, like radius).

Change

::

    PA = params[0 + offsetIndex];
    ell = params[1 + offsetIndex];
    I_0 = params[2 + offsetIndex];
    sigma = params[3 + offsetIndex];

to

::

    PA = params[0 + offsetIndex];
    ell = params[1 + offsetIndex];
    I_0 = params[2 + offsetIndex];
    fwhm = params[3 + offsetIndex];
    beta = params[4 + offsetIndex];

Then, at the end of the method, replaced this line

::

    twosigma_squared = 2.0 * sigma*sigma;

with this (which computes the "alpha" parameter of the Moffat function)

::

    double exponent = pow(2.0, 1.0/beta);
    alpha = 0.5*fwhm/sqrt(exponent = 1.0);

D. Changes to the ``CalculateIntensity`` method:

XXX

Other Potential Issues
----------------------

If your new image function has an analytic expression for the total
flux, then you might consider overriding the CanCalculateTotalFlux
method to return ``true`` and then override the ``TotalFlux`` method so
that it calculates and returns the total flux. (The default is to let
``makeimage`` estimate the total flux numerically, by generating a large
image using the image function and summing all the pixel values.)
