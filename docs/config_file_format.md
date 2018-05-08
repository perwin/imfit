# Configuration File Format

Imfit and makeimage *always* require a configuration file, a text file which describes the
model to be fit to the data (in the case of imfit) or to be generated (in the
case of makeimage), plus some optional information about the image itself.

This page provides a basic description of the configuration file format.

### Blank lines and Comments

Any blank lines in the file are ignored. Any lines beginning with the hash/pound
symbol "#" are also ignored.

Finally, if a hash symbol is found in the middle of a line, then it and the
rest of the line are ignored.  This allows comments to be added on the same
line as a meaningful entry (function declaration, parameter name and value, etc.).


## (Optional) Prelude: Describing the Image

The configuration file can start with a prelude which provides
information about the data image (for imfit) or the output model image (in the
case of makeimage). This is in the form of lines containing single "NAME value" pairs, e.g.   

    GAIN 4.5

These do not *have* to be provided, and you can provide as many or as few of
these as are needed (the exception being the NCOLS and NROWS entries for a
makeimage configuration file -- both of those should be provided).

The allowed entries for an imfit configuration file are:

* GAIN -- the A/D gain value (electrons/count) for the image

* READNOISE -- any Gaussian read noise (electrons)

* EXPTIME -- total integration time *if* pixel values are counts/sec

* NCOMBINED -- number of images combined *if* pixel values are mean or median

* ORIGINAL_SKY -- any original constant sky-background value which was subtracted from the image

The allowed entries for a makeimage configuration file are:

* NCOLS -- integer defining the width (number of columns) of the output image

* NROWS -- integer defining the height (number of rows) of the output image

All of this information can *also* be provided via command-line options, which if present
will override any corresponding values in the configuration file.


## Main Section: Defining the Model

The main (and only required) section of the configuration file is the part which
describes the actual model to be fit (or generated and saved in the case of makeimage).
This consists of one or more **function blocks**.

### Function Blocks

A **function block** is a central-coordinate (x0,y0) specification followed by one or
more **image-function declarations**. This specifies a set of image functions which
share the same center.

    X0  <initial_value>   [<limits>]
    Y0  <initial_value>   [<limits>]
    FUNCTION <function_name>
    <parameter specifications ...>

The optional limit specifications for X0 and Y0, which are only used by imfit, can have one of two forms:
1. The word "fixed", which indicates that the parameter should be held fixed
2. A comma-separated pair of numbers specifying lower and upper limits

### Image-Function Declarations

An **image-function declaration** consists of the word "FUNCTION" followed by the
name of an image function. (Recall that you can use `--list-functions` to
get a list of the available image functions and `--list-parameters` to get that
list along with the parameter list for each function.) This is followed by the
list of parameter values; these are the values for generating a model image with
makeimage, or the *initial* values for the fitting process
in imfit.

    FUNCTION <function_name>
    <parameter_name1>    <initial_value>   [<limits>]
    <parameter_name2>    <initial_value>   [<limits>]
    [etc.]

A simple example:

    FUNCTION Exponential
    PA    90
    ell   0.5
    I_0   100
    h     15

The optional limit specification for the parameters are exactly like those for the
X0 and Y0 values (see above):
1. The word "fixed", which indicates that the parameter should be held fixed
2. A comma-separated pair of numbers specifying lower and upper limits

Note that limit specifications are actually *required* by the Differential Evolution
solver; they are optional for the other solvers (including the default Levenberg-Marquardt solver).
Limit specifications are *ignored* by makeimage.


## Examples

### Example: Single-function image

This is based on configuration file in the `examples/` subdirectory. If
used with imfit, it fits a single Sersic function to a 256 by 256-pixel
cutout of an SDSS image. If used with makeimage, it generates a 2D Sersic-function
image using the specified parameters. This example does *not* contain gain or read-noise
specifications (for use with imfit). Lower and upper limites are listed for all parameters
for use with imfit.

    X0   129.0    125,135
    Y0   129.0    125,135
    FUNCTION Sersic
    PA    18.0    0,90
    ell    0.2    0,1
    n      1.5    0,5
    I_e    15     0,500
    r_e    25     0,100


### Example: Single function block with two functions

This is a modification of the previous configuration file, using an Exponential function
along with the Sersic function. In addition, the Sersic index *n* is held fixed with a
value of 4 (making the Sersic profile a de Vaucouleurs profile). Both functions share the 
same center, and are thus part of a single function block. This version also includes
an image-description prelude.

    GAIN          4.725
    READNOISE     4.3
    ORIGINAL_SKY  130.14
    
    X0   129.0    125,135
    Y0   129.0    125,135
    FUNCTION Sersic
    PA    18.0    0,90
    ell    0.2    0,1
    n      4      fixed
    I_e    15     0,500
    r_e    25     0,100
    FUNCTION Exponential
    PA    18.0    0,90
    ell   0.5     0,0.8
    I_0   100     1,500
    h     50      5,500

### Example: Multiple function blocks

TBD.

