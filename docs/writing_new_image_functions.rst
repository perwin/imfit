How to Write New Image Functions
================================

Overview
--------

Imfit comes with a variety of 2D surface-brightness functions (“image
functions”). But it is designed to make it relatively simple to add
*new* image functions, to describe objects or substructures which are
not well modeled by the existing functions. This is done by writing
additional C++ code and re-compiling the various Imfit programs
(``imfit``, ``imfit-mcmc``, ``makeimage``) to include the new function
or functions.

This page provides some guidelines on how to write new image functions
(see also the “Rolling Your Own Image Functions” section of the `Imfit
manual <https://www.mpe.mpg.de/~erwin/resources/imfit/imfit_howto.pdf>`__,
from which this page is excerpted). Even if you’re not very familiar
with C++, it should hopefully not be too difficult to write a new image
function, since it can be done by copying and modifying one of the
existing image functions.

A new image function is written in C++ as a subclass of the
FunctionObject base class, which is declared and defined in the
``function_object.h`` and ``function_object.cpp`` source-code files (in
the ``function_objects/`` subdirectory of the source-code distribution).

A new image-function class should provide its own implementation of
(i.e., “override”) the following public FunctionObject methods, which
are defined as virtual methods in the base class:

-  The class constructor – in most cases, the code for this can be
   copied from any of the existing FunctionObject subclasses, unless
   some special extra initialization is needed.

-  ``Setup()`` – this is used by the calling program to supply the
   current set of function parameters (including the (*x*\ 0,\ *y*\ 0)
   pixel values for the center of the function set) prior to determining
   intensity values for individual pixels. The parameter values in the
   input array should be assigned to the appropriate data members of the
   class. This is also a convenient place to do any general calculations
   which depend on the parameter values but don’t depend on the exact
   pixel (*x*,\ *y*) values.

-  ``GetValue()`` – this is used by the calling program to obtain the
   surface brightness for a given pixel location (*x*,\ *y*). In
   existing FunctionObject subclasses, this method often calls other
   (private) methods to handle details of the calculation.

-  ``GetClassShortName()`` – this is a class function which returns the
   short version of the class name as a string.

The new class should also redefine the following internal class
constants:

-  ``N_PARAMS`` — the number of input parameters (*excluding* the
   central pixel coordinates);
-  ``PARAM_LABELS`` — a vector of string labels for the input
   parameters;
-  ``FUNCTION_NAME`` — a short string describing the function;
-  ``className`` — a string (no spaces allowed) giving the official name
   of the function.

The ``add_functions.cpp`` file should then be updated by:

-  including the header file for the new class;
-  adding 2 lines to the ``PopulateFactoryMap()`` function to add the
   ability to create an instance of the new class. (Look for the comment
   line “// ADD CODE FOR NEW FUNCTIONS HERE” in the file.)

Finally, the name of the C++ implementation file for the new class
should be added to the ``SConstruct`` file to ensure it gets included in
the compilation; the easiest thing is to add the file’s name (without
the ``.cpp`` suffix) to the multi-line ``functionobject_obj_string``
string definition. (Or look for the comment line “# ADD CODE FOR NEW
FUNCTIONS HERE” in the ``SConstruct`` file.)

Existing examples of FunctionObject subclasses can be found in the
``function_objects/`` subdirectory of the source-code distribution, and
are the best place to look in order to get a better sense of how to
implement new FunctionObject subclasses.

A Simple Example
----------------

To demonstrate the basics of writing a new image function, we’ll modify
the existing Gaussian class to make a class called NewMoffat, which
produces an elliptical structure with a Moffat radial profile. (Such a
function already exists in Imfit under the name Moffat, so this is
really a redundant exercise.)

Three basic changes to the existing ``func_gauss.h/cpp`` files are
needed:

1. Change the class name (in this case, from “Gaussian” to “NewMoffat”);

2. Change the code which actually computes the function;

3. Add, rename, and delete class data members to accommodate the new
   algorithm.

1. Create and Edit the Header File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming we’re in the Imfit source code directory, ``cd`` to the
``function_objects/`` subdirectory, copy the header file
``func_gaussian.h`` to ``func_new-moffat.h``, and edit the file to
change the following lines:

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

   const int N_PARAMS = 4;
   const char PARAM_LABELS[][20] = {"PA", "ell", "I_0", "sigma"};
   const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "counts/pixel", "pixels"};
   const char FUNCTION_NAME[] = "Gaussian function";

to

::

   #include "func_new-moffat.h"

   const int N_PARAMS = 5;
   const char PARAM_LABELS[][20] = {"PA", "ell", "I_0", "fwhm", "beta"};
   const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "", "counts/pixel", 
                               "pixels", ""};
   const char FUNCTION_NAME[] = "Moffat function";

B. In the remainder of the file, change all references to the class name
from ``Gaussian`` to ``NewMoffat`` (e.g., ``Gaussian::Setup`` becomes
``NewMoffat::Setup``).

C. Change the ``Setup`` method. Here, you’ll need to change how the
input parameter array is converted into individual parameters, and do
any useful pre-computations (i.e., computations that depend on the
parameter values, but not on individual pixel values or values derived
from the latter, like radius).

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

with this (which computes the “alpha” parameter of the Moffat function)

::

   double exponent = pow(2.0, 1.0/beta);
   alpha = 0.5*fwhm/sqrt(exponent = 1.0);

D. Changes to the ``CalculateIntensity`` method:

Although it is the public method GetValue which is called by other parts
of the program, we don’t actually need to change the current version of
that method in this example. The code in the original Gaussian version
of GetValue converts pixel positions to a scaled radius value, given
input values for the center, ellipticity, and position angle, and then
calls the private method CalculateIntensity to determine the intensity
as a function of the radius. Since we’re still assuming a perfectly
elliptical shape, we can keep the existing code. (GetValue also includes
possible pixel subsampling, which is useful for cases where intensity
changes rapidly one scales of a single pixel; we’ll apply a simple
modification for the Moffat function later on.)

So in this case we actually implement the details of the new function’s
algorithm in CalculateIntensity. Replace the original version of that
method with the following:

::

   double NewMoffat::CalculateIntensity( double r )
   {
     double  scaledR, denominator;

     scaledR = r / alpha;
     denominator = pow((1.0 + scaledR*scaledR), beta);
     return (I_0 / denominator);
   }

E. Changes to the ``CalculateSubsamples`` method:

Although pixel subsampling is performed in the GetValues method, the
determination of whether or not to actually \*do*\* the subsampling –
and how much of it to do – is determined in CalcualteSubsamples.

For the Gaussian function, subsampling can be useful happen when *r* < 1
*and* sigma < 1. The equivalent for the Moffat function would be *r* < 1
and alpha < 1, so change the line in CalculateSubsamples that says

::

   if ((sigma <= 1.0) && (r <= 1.0))

to say

::

   if ((alpha <= 1.0) && (r <= 1.0))

At this point, most of the work is done. We only need to update the code
in ``add_functions.cpp`` so it knows about the new function and update
the ``SConstruct`` file so that the new function is included in the
compilation.

Adding Parameter Unit Definitions
---------------------------------

Imfit version 1.9 added the possibility of including strings describing
the intended units for certain parameters of the image functions. This
is purely cosmetic – the image functions ignore them – but it can be
useful for the user as a reminder of what individual parameters mean.
(Example: the Exponential image function has the units label “deg (CCW
from +y axis)” for the PA (position angle) parameter, “counts/pixel” for
I_0, and “pixels” for the scale length parameter h.)

You do not *have* to provide this information in any image functions you
write. If you wish to, you need to properly define the ``PARAM_UNITS``
constant at the top of the source file (see 2.A above for an example);
assign the unit-description strings (the empty string – “” – is used for
an parameters without units) to the ``parameterUnits`` vector in the
constructor (see any of the standard image functions for an example of
how to do this), and set the ``parameterUnitsExist`` data member to
``true`` (again, see the constructor of any of the standard image
functions for examples).

Other Potential Issues
----------------------

If your new image function has an analytic expression for the total
flux, then you might consider overriding the ``CanCalculateTotalFlux``
method to return ``true`` and then overriding the ``TotalFlux`` method
so that it calculates and returns the total flux. (The default is to let
``makeimage`` estimate the total flux numerically, by generating a large
image using the image function and summing all the pixel values.)

If your new image function is meant to represent the image *background*
(as in the case of the built-in functions FlatSky and TiltedSkyPlane),
then you may not want ``makeimage`` trying to calculate the “total flux”
for the component. In this case, you should add a line to your class’s
constructor which sets the data member ``isBackground`` to ``true`` (as
is done in, e.g., ``func_flatsky.cpp``).
