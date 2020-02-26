# Imfit Frequently Asked Questions

This is preliminary version of a page for assorted (potentially) Frequently
Asked Questions about using Imfit.


## How do I include a true point source (star, AGN, etc.) in a model?

Although point sources can be approximated by something like a Gaussian function
with a very small size (e.g, 0.1 pixels), the best way to include point sources
is with the PointSource function:

    FUNCTION PointSource
    I_0		<initial value  [limits]>

This uses the user-supplied PSF image, shifting it (via bicubic interpolation) and
scaling its brightness using the `I_0` parameter (which corresponds to the
total flux of the source).


## How can I determine B/T (bulge/total) and other flux ratios?

Assuming you have *either* an `imfit`/`makeimage` config file *or* the
best-fit-parameters output file produced by running `imfit` (which has
the same general format), you can run `makeimage` in "--print-fluxes"
mode, which will compute the fluxes of the individual components and
their fractions of the total flux of the model:

    $ makeimage config_or_bestfit_file.dat --print-fluxes [--nosave]

(The "--nosave" option tells `makeimage` to skip saving the model image;
this is useful if you just want to know what the fluxes and ratios are.)

The output will look something like this:

    Component                 Flux        Magnitude  Fraction
    Sersic                    1.0800e+06      ---     0.09038
    Exponential               4.8448e+06      ---     0.40545
    Sersic_GenEllipse         5.6323e+05      ---     0.04713
    Exponential               5.4612e+06      ---     0.45703

    Total                     1.1949e+07

Note that the flux units are in counts (i.e., the units of individual
pixel values). To get outputs in magnitudes as well, see the next question.

By default, `makeimage` computes component fluxes by generating a 5000 x
5000 pixel internal image and summing up the pixel values (a few
functions -- PointSource, Gaussian, Exponential, and Sersic -- are able to
compute their total fluxes analytically, which is faster). Should you be
working with images larger than this size -- or with models whose flux will
extend significantly outside the default size -- you can specify a larger
(or smaller) summation image size with the "--estimation-size=<int_value>",
where the specified value is the size of the (square) image.

Finally, fluxes for "background" image functions (e.g., FlatSky), which
are unbounded and would just scale with the summation image size, are
treated as zero.


## How can I determine magnitudes of model components?

Imfit works in units of counts/pixel, so it does not by default worry
about magnitude systems, zero points, etc. However, you can ask the
`makeimage` program, as part of the "--print-fluxes" mode, to compute
magnitudes for individual components (and the whole model):

    $ makeimage config_file.dat --print-fluxes --zero-point=<ZP>

The "ZP" value will be used to compute the magnitude of each component, as

    m = ZP - 2.5 log10(flux)

The output will look something like this (using "--zero-point=25"):

    Component                 Flux        Magnitude  Fraction
    Sersic                    1.0800e+06   9.9164     0.09038
    Exponential               4.8448e+06   8.2868     0.40545
    Sersic_GenEllipse         5.6323e+05   10.6233     0.04713
    Exponential               5.4612e+06   8.1568     0.45703

    Total                     1.1949e+07   7.3066
