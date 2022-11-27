# Change Log for Imfit

(Formatting and design based on Olivier Lacan's [Keep a CHANGELOG](http://keepachangelog.com/))

## 1.9.0 -- 2022-11-27
### Added:

- Imfit will attempt to identify how many *physical* CPU cores your
computer has, and then set the number of threads equal to that. In some
cases (particularly with smaller images), this can make PSF convolution
faster than the previous default of setting the number of threads equal
to the number of *logical* cores, which usually includes so-called
symmetric multithreading (aka "hyperthreading") and turns out to (sometimes)
make PSF convolutions slower.

- New image function: PointSourceRot. This is identical to PointSource,
except with an extra "PA" parameter which specifies CCW rotation of the
PSF image. Thanks to Alex Borlaff for suggesting this.

- Unit strings for image-function parameters: the standard image functions
will now print unit strings (e.g., "pixels", "counts/pixel", "deg")
after the parameter names (when imfit or makeimage is called with
"--list-parameters") and after the values (in imfit's output). The
printing of parameters without units (e.g., ellipticity, Sersic n) is 
unchanged. (You can optionally add unit strings to the printout for
your own image functions; this is explained in Section 14 of the 
Imfit manual.)

- The PointSource (and PointSourceRot) image functions now include the
option of Lanczos2 instead of the default bicubic interpolation. This
is enabled by adding the following three lines immediately after the 
"FUNCTION PointSource" line in a configuration file: "OPTIONAL_PARAMS_START",
"method   lanczos2", "OPTIONAL_PARAMS_END".

- It is now possible to write image functions which, like the updated
versions of PointSource and PointSourceRot, take one (or more)
*optional* parameters to specialize the function. These are *not* like
the standard parameters -- which are changed by the optimization
algorithm as part of the fitting process -- instead, they are only
used/applied at startup. These come in the form of key-value pairs (the
key is the name of the parameter), with both specified as text strings
in the configuration file. If supplied, they are passed to the image
function after it is initialized but before any fitting; how they values
are used is up to the individual function (e.g., they can be converted
to floating-point values, used to select one of several internal
algorithms, interpreted as the filename of an extra configuration file
with options and/or data for the function, etc.).
    

### Changed:

- If you type Control-C to interrupt a fit (you think it's taking too long,
you decide you specified the model wrong, you forgot to specify a PSF image, etc.),
`imfit` now halts, prints the *current* best parameter values, *and*
saves them in a special output file ("current_parameters_imfit.dat").

- The `--print-fluxes` option for makeimage now automatically turns off
image saving, so you no longer have to *also* specify `--nosave`. (This
is because, in my experience at least, one almost always uses
`--print-fluxes` just to see the component fluxes and fractions, and so
it's tedious to have to remember to also specify `--nosave`.)



## 1.8.0 -- 2020-07-26
### Added:

- Function "labels" in configuration files: you can now annotate individual
functions in config files (e.g., "FUNCTION Sersic # LABEL bulge");
these labels will be propagated to output best-fit parameter files and
`--print-fluxes` and `--save-fluxes` output. Everything after the keyword
"LABEL" (must be in all-caps) is the "label"; note that the comment character "#"
is required!

- New image function: TiltedSkyPlane. This represents a linear 2D
gradient background (a "tilted plane") as the image background. Based on
the InclinedFlatSky function suggested by Dan Prole (danjampro).

- New image function: GaussianRingAz. This is similar to GaussianRing
(an elliptical ring with a Gaussian radial intensity profile), but has a
peak intensity that varies with azimuth around the ring (see Erwin et
al. 2021 for its use in modeling images of the barred galaxy NGC 4608).

- New image function: FlatBar. This is an ad-hoc representation of the
outer, vertically thin part of a two-component, vertically buckled bar
in a massive galaxy, seen at low to moderate inclinations. It is meant to be
used in conjunction with, e.g., a Sersic or Sersic_GenEllipse component (which
would represent the inner, vertically thickened "boxy/peanut-shaped bulge"
part of the bar. (See Erwin et al. 2021 -- especially Appendix~A -- for its use
in modeling images of two barred galaxies).

- New option for makeimage: `--save-fluxes`, which saves the
`--print-fluxes` component-fluxes-and-fractions output of makeimage to a
text file. (Thanks to Iskren Georgiev for requesting this.)

- New option for fitting using Differential Evolution: `--de-lhs`, which
uses Latin hypercube sampling to set up the initial set of trial
parameter vectors, instead of the default pure-uniform sampling. (This
may in some cases more evenly sample the parameter volume, though
preliminary tests suggest it may be no better than the default.)

### Changed:

- The documentation now (mostly) refers to "function blocks" (i.e., one
or more image functions in a configuration file that share the same
central pixel coordinates X0,Y0) as *function sets*, to match the usage
in PyImfit better.

- Compilation of the Imfit programs (`imfit`, `imfit-mcmc`, and `makeimage`)
now *requires* the GNU Scientific Library; i.e., it is no longer "optional
but strongly recommended".





## 1.7.1 -- 2019-08-07
### Changed:

- The imfit_funcs.py module now calculates the Sersic b_n parameter correctly
for all values of n > 0, and a dependence on the mpmath module has been
removed.

### Fixed:

- An error in the SConstruct file introduced in 1.7.0 meant that OpenMP
multi-threading was effectively turned off, making Imfit slow. This
affected all pre-compiled versions *and* anything compiled from source
with the SConstruct file. This has been fixed -- all users of v1.7.0
are encouraged to update to restore Imfit's traditional speed!



## 1.7.0 -- 2019-07-10
### Added:

- New image function: FerrersBar2D. This is an elliptical 2D image
function (with boxy/disky generalized-ellipse isophotes) based on the
Ferrers ellipsoid, but with the latter's 3D density value interpreted as 2D
surface-brightness. (The true 3D Ferrers ellipsoid has been available as
FerrersBar3D since version 1.5.)

- You can now compile a static library file (`libimfit.a`) that
potentially allows using most of Imfit's functionality from other
programs (this does not include the MCMC code, nor does it include the
CFITSIO-based reading and writing of FITS files). This is primarily to
support PyImfit (the Python wrapper for Imfit); there isn't any API
documentation to speak of, so it's not really useful for anything else
at this point.

### Changed:

- When doing bootstrap resampling, a one-line progress bar (with updating
iteration count and percentage completion) is drawn to stdout; this replaces
the printing of individual iteration numbers in previous versions, which
could end up spilling over many lines of output. (Thanks to Iskren Georgiev
for suggesting this.)

- The SConstruct file has been updated to allow compiling on a Mac using
Apple's Clang compiler, which (apparently since late 2017) finally has
support for OpenMP! This *does*, however, require installing`libomp` (e.g., `brew
install libomp`) and using the `--clang-openmp` flag when compiling; the
documentation has been updated to describe this.

- The precompiled macOS binary is now built with Apple's Clang compiler,
rather than GCC. (This makes the binaries somewhat smaller, since they don't
have to include statically linked copies of the GCC libraries.)

### Fixed:

- The PointSource function was not producing the correct results *if* it was
located within an oversampling region. This has been fixed. (Thanks to Iskren
Georgiev for spotting this.)

- The output of bootstrap-resampling for imfit for fits with image sections would
not include the correct (full-image) X0,Y0 values at the beginning of each parameter
line. Note that this was only a problem for the printed-to-screen information;
saving the correct X0,Y0 values to the bootstrap-output file was fixed back in
1.5.0. (Thanks to Chris Pritchett for spotting this.)

- PSF images whose pixels values sum to <= 0 are now caught as input errors
(in the default case of automatic PSF-image normalization; if the `--no-normalize`
flag is used, then no checks of the PSF image other than for NaN values are done).



## 1.6.1 -- 2018-11-15
### Changed:

- The output-file-reading code in imfit.py will now use pandas.read_csv instead of
numpy.loadtxt, if pandas is installed. This is significantly faster when reading in
large MCMC output files. (Thanks to Justus Neumann and Iskren Georgiev for
suggesting this.)

- Library code included in the compiled binaries has been updated to version 2.5 of
the GNU Scientific Library, version 3.3.8 of FFTW, version 3.450 of cfitsio, and 
version 2.5.0 of NLopt.

### Fixed:

- imfit-mcmc would use the wrong parameter limits for X0 and Y0
parameters if an image section was specified (e.g.,
"<data-file>.fits[x1:x2,y1:2]"); in some cases this could result in
failed convergence or mis-constrained results. (Thanks to Iskren Georgiev for
spotting the problem.)



## 1.6.0 -- 2018-03-24
### Added:

- New image function: PointSource. This is an interpolating function
which takes the user-supplied PSF image, rescales it in intensity, and
interpolates it to the position specified by the X0 and Y0 parameters of
a function block. (If an oversampled PSF is supplied, then that will be
used for the specified oversampled region or regions.)
This is more accurate for true point sources than the current approach
of constructing a narrow Gaussian and convolving that with the PSF
image. The interpolation is done with the GNU Scientific Library bicubic interpolation
function; future version may offer the option of something like
interpolation with a Lanczos kernel as well. (Thanks to Corentin Schreiber for
suggesting this development and helping test it!)

### Fixed:

- Convolution with oversampled PSFs was not being done when single-function
images were generated (i.e., by makeimage with the --output-functions option).
This has been fixed.

- Previously, parameters in config file for MCMC mode (imfit-mcmc) which were
missing lower and upper limits would silently convert to fixed values of 0.
Now, when parameters without limits are supplied, an error is signaled and
imfit-mcmc quits. (Thanks to Corentin Schreiber for reporting this.)

- The bash tab-completion file introduced in 1.5.0 did not include filename
completion (i.e., if you typed "imfit " and then the beginning of a filename
in the current directory and *then* pressed the TAB key, the bash shell's
normal filename completion wouldn't work). This has been fixed.

- Pre-compiled versions of the code now link to version 3.430 of the CFITSIO
library, which includes important security fixes.


## 1.5.0 -- 2017-08-16
### Added:

- Multiple oversampling regions: Imfit can now optionally convolve
multiple subsections of the model image with oversampled PSFs, instead
of just one subsection as in previous versions. (This applies to all the
individual programs: imfit, imfit-mcmc, and makeimage.)<br><br>
By adding multiple invocations of `--overpsf_region x1:x2,y1:y2`, each
region will be convolved with the same oversampled PSF (specified by
`--overpsf <fits-filename>` and `overpsf-scale <int>`). Alternately, you
can also add multiple invocations of `--overpsf <fits-filename>` and
`overpsf-scale <int>` to specify that each oversampled region should be
convolved with its own, separate PSF image.

- New image function: FerrersBar3D. This is the classic Ferrers (1877)
ellipsoid often used in modeling of orbits in barred-galaxy potentials,
which implements a bounded triaxial ellipsoid seen at arbitrary orientation
with luminosity density following the Ferrers mass-density
form. As with the ExponentialDisk3D and GaussianRing3D functions, this
uses integration along the line of sight to construct the model image.

- Command-line flag `--no-normalize` will tell imfit, imfit-mcmc, and makeimage
to *skip* normalization of any input PSF images (standard or oversampled). 
(Thanks to Corentin Schreiber for requesting this.)

- Tab auto-completions (for the Bash shell)! Imfit now includes a shell
script (`extras/imfit_completions.bash`) which can be used to define
auto-completions of imfit/imfit-mcmc/makeimage command-line options with
the TAB key when using the Bash shell.

### Changed:

- Output MCMC chains produced by `imfit-mcmc` are now numbered starting with 1
instead of with 0. (This is more consistent with the rest of Imfit, which uses
1-based numbering almost everywhere outside of the actual C++ code. My apologies if
this messes up anyone's analysis code!)

- I have started (very slowly) updating some of the code to make use of
C++-11 features. This should have no visible effects, but might start to
be relevant for people hacking around with the code. It also means that
Imfit now requires compilers that support C++-11 (e.g., GCC versions 4.8
or newer).

- The "--nosubsampling" flag has been renamed to "--no-subsampling", to be
more consistent with other multi-syllabic option names.

- Memory-use estimation now includes the effects of multiple PSF oversampling
regions.

- Binaries for MacOS X (or "macOS") are now compiled with GCC 7 instead of GCC 5.
(But if you're compiling the source code yourself, any version of GCC from 4.8
onward will do.)

- I am no longer including pre-compiled versions for Mac OS X 10.6 or 10.7.

### Fixed:

- When bootstrap-resampling results (imfit) or MCMC chains (imfit-mcmc)
were saved in output files, the X0,Y0 coordinates were erroneously *not*
corrected if an image section was specified (i.e., fitting
`data.fits[x1:x2,y1:y2]` as opposed to `data.fits`); X0 and Y0 values
are now corrected back to full-image coordinates (as was already done
for printed and saved best-fit parameter values). (Thanks to Iskren
Georgiev for spotting this in the MCMC output.)

- Printouts and output files listing the best-fitting parameter values
were supposed to include a blank line between function blocks (for easier
reading), but this only happened when Levenberg-Marquardt fits were
done. This now works for fits with N-M simplex and Differential Evolution
as well.

- In describing the PSF oversampling options, the documentation
(erroneously) gave an example of `--overpsf_region [x1:x2,y1:y2]`, when
the correct format was `--overpsf_region x1:x2,y1:y2`. The documentation
has been updated to recommend the latter format; however, the code has
also been changed so that *both* formats are now acceptable.
(Thanks to Iskren Georgiev for reporting this.)

- The binary Mac distribution of version 1.4 was accidentally linked with
a version of the FFTW library which did *not* include SSE2 vector math
support; consequently, fits with PSF convolution were about 5-10% *slower*
than in version 1.3. The binary distribution for Macs is now linked to
an FFTW library with both SSE2 and AVX support, and so should be as
fast as version 1.3 again (and possibly a little faster on recent CPUs).

- The Python function python/imfit.py was missing an "import numpy as np"
statement; this has been fixed.


## 1.4.0 -- 2016-11-10
### Added:

- Imfit (the package) now includes a new program for doing Markov chain
Monte Carlo (MCMC) analysis of Imfit models: "imfit-mcmc". Rather than
fitting the model to the data, this estimates the posterior-probability
(i.e., the likelihood, equivalent to chi^2 in the default case)
distribution for the model given the data. This uses the DREAM
(DiffeRential Evolution Adaptive Metropolis) algorithm of Vrugt et al.,
an efficient multiple-ensemble approach. See Chapter 13 of the
documentation for more details.

- There are now some very simple Python functions in python/imfit.py which can read
in output from bootstrap resampling (imfit.GetBootstrapOutput) and from MCMC chains
(imfit.GetSingleChain, imfit.MergeChains).

- Makeimage now has a new option "`--timing <N>`", which tells the program to
generate the model image *N* times and compute the mean time taken per image. 
(No images are saved.) This is potentially useful if you want an approximate
idea of how long a fit (or MCMC run) might take.

- Imfit now has a new option "`--seed <N>`", which allows the user to
specify a particular seed (positive integer) for the random number
generator. This is mainly useful if you want to test the Differential
Evolution solver, bootstrap resampling, or the MCMC program by forcing
them to use the same sequence of pseudo-random numbers. (The default
behavior is to always initialize the RNG with the current system time.)


### Changed:

- Tweaks to the Differential Evolution coefficients now make DE fits about
30--80% faster.

- Imfit's internal PSF normalization (which takes place when PSF images are
read in during startup) now uses the Kahan summation algorithm, which gives
more accurate results when many pixel values are very close to zero. This should
generally have no effect at all, but might be relevant when using very large
PSF images with very faint, extended wings.

- Some of the image functions now have the ability to calculate the
total flux analytically rather than by the standard "integrate over
large internal image" method. This makes estimating total fluxes via
makeimage's "`--print-fluxes`" option much faster for those functions. In
some cases, the resulting fluxes are different from the previous
approach, but only at levels <~ 10^-4 (i.e., about 0.0001 mag). Affected
functions: Gaussian, Exponential, Sersic (the latter only if
makeimage was compiled with the GNU Scientific Library).

- Memory-use estimation now prints values < 1 MB as KB. (Which is pretty
unnecessary, except that "0.0 MB" looked kind of silly as an output when
the estimated memory use was < 100 KB.) Memory estimation is also slightly
more accurate now if one or more parameters are held fixed and L-M minimization
is being used.

- The SConstruct file no longer checks to see if there is a version of
Xcode earlier than 5 installed on a Mac (in which case it would try to
use llvm-g++-4.2 as the compiler), since that's increasingly
ancient. The choices for compiling on a Mac are now basically: 1) Use
GCC installed via a package manager (e.g., Homebrew, Fink, MacPorts); or 2) Use Apple's version of clang without OpenMP support.



## 1.3.1 -- 2016-05-05
### Added:

- The full codebase for Imfit is now available on Github at
[github.com/perwin/imfit](https://github.com/perwin/imfit), with
automated unit and build tests executed via Travis CI. (The compact,
all-you-need-to-compile version of the source code will continue
to be available via the main web page.)

### Changed:
- Any pixels in the noise/error image which have non-finite values
(NaN, +/- infinity) are now automatically masked, rather than causing
Imfit to quit with an error message. (Thanks to Dave Wilman for requesting this.)

- Any pixels in the *mask* image which have non-finite values (NaN, +/- infinity) are now
treated as part of the mask; previously, they were (unintentionally) treated as
indicating *valid* data pixels. (This is arguably more of a "fix" than a "change",
since it seems unlikely anyone would intend NaN mask values to actually indicate
good data.)

- FITS headers for output images are slightly more informative (e.g.,
headers for best-fit model image and residual image now include name of original
data image.

- I have dropped the pre-compiled 32-bit Mac OS X version from the standard
binary distributions. (If you have a 32-bit-only machine running OS X 10.6 or
10.7, it should still be possible to compile from source code.)


### Fixed:
- Better checking of input images: FITS files which do *not* have a valid 2D image
in their primary header-data unit are now recognized and rejected. (Previously,
they were sometimes rejected with a rather opaque error message and sometimes read in;
the latter would, obviously, cause problems....)

- Added checking for memory-allocation failures.

- Fixed minor bug which could cause segfault at very end of imfit's run
*if* bootstrapping was done without requesting that bootstrap output be
saved (but only on Linux).



## 1.3 -- 2015-10-05
### Added: 
- New image function: modified (or "empirical") King profile, as
described in Elson et al. (1999) and Peng et al. (2010). This comes in
two versions: one which specifies the tidal/truncation radius as a free
parameter (ModifiedKing), and the other which has the concentration as a
free parameter in place of the tidal radius (ModifiedKing2).

- Imfit and makeimage can now export sample configuration files, using the "`--sample-config`"
command-line flag.

- The parameter uncertainties generated by the Levenberg-Marquardt minimization algorithm 
are now recorded in the output best-fit-parameters file in the same way that they are 
printed to the screen; e.g., "X0              32.9439 # +/- 0.0128".
(Thanks to Rebecca Lange and Semyeong Oh for suggesting this.)

- Imfit now reports the time taken at the end of its run (including separate times for
fitting and bootstrap resampling), unless the "`--silent`" command-line flag is used.

- The Imfit web pages now include a simple tutorial.

### Changed:
- Import of external error/noise/weight-map images has been changed, but only in the case of
"`--errors-are-weights`" option; previously, these were treated as though they were
identical to imfit's internal 1/sigma weights, but *described* in the documentation
and paper as though they were 1/sigma^2. If error map *is* specified as *weights*,
it is now read in assuming it has 1/sigma^2 values, and saved weight maps (via the
"`--save-weights`" option) are converted to the same format. (Thanks to Semyeong Oh for
spotting this problem.)

- PSF convolution code has been rewritten to use smaller arrays and
slightly faster algorithms (taking advantage of the specialized
real-to-complex, complex-to-real transformations in the FFTW library).
The result is a reduction in memory use by ~ 20%, and a speedup in doing
PSF convolutions (see below).

- Precompiled binaries now use a version of the FFTW library
which includes SSE2 vectorization (available in all Intel and AMD x86-type
CPUs manufactured since about 2003).

- The combination of the preceding two changes appears to make fits with
PSF convolution about 30% faster on average than they were in version
1.2 (in some cases they can even be twice as fast).

- Printouts and output files listing the best-fitting parameter values (and errors) 
now include a blank line between function blocks, for easier reading.

- Output from "`imfit --help`" has been reorganized into a (hopefully) more logical form.

- Memory-use estimates now account for convolutions with oversampled PSFs (and the
changes in the convolution code).

- Updates to documentation.

### Fixed:
- Convolutions when the model contains negative pixel values are now
done correctly. (Thanks to Semyeong Oh for spotting the problem.)



## 1.2.1 -- 2015-06-30
### Fixed:
- Fixes a bug in the bootstrap-resampling summary output (if none of the parameters
had user-specified limits, then "[fixed]" would erroneously be printed in place of 
the actual bootstrap confidence intervals). (Thanks to Semyeong Oh for spotting this.)

- Fixes compilation bug with "`--no-nlopt`" option. (Thanks to Semyeong Oh for spotting this.)

- Added note to imfit_howto.pdf explaining how weights are internally handled, which
disagrees with how they are presented in the paper (this is largely cosmetic, except
when using the "`--errors-are-weights`" and "`--save-weights`" option, and will be resolved in v1.3).



## 1.2 -- 2015-05-27
### Added: 
- Imfit can now optionally convolve part of the model image with an
oversampled PSF. E.g., you can specify that a 10x10-pixel region
centered on a galaxy nucleus should be modeled using a
five-times-smaller pixel size and convolved with a corresponding
(five-times-oversampled) PSF image. The resulting oversampled and
convolved sub-image is then downsampled back to the main image pixel
scale before the model is compared with the data image (This is in
addition to the standard PSF convolution that imfit already allows,
where the PSF image and the computed model image have the same pixel
scale as the data image.)

- The saved best-fit parameter file produced by imfit now include a
brief summary of the fitting process and its outcome (fitting statistic
and minimization algorithm used, final best-fit value of statistic, AIC,
BIC), written in the header of the file as comments. (Thanks to Colleen Gilhuly
for suggesting this.)

- Imfit now attempts to estimate (and print) the total amount of memory
needed before it starts the fitting process. This estimate is crude (and
purely advisory), but may be useful in cases when fitting very large
images (and/or using very large PSF images) might run up against your
computer's memory limits.  Note that this currently does *not* account
for any use of an oversampled PSF. (Thanks to Lee Kelvin for helping demonstrate
the utility of this.)

### Changed:
- Updates to documentation.


### Fixed:
- Better error detection for mangled parameter lines in config files (which could
cause mysterious segmentation faults otherwise).



## 1.1 -- 2014-10-22
### Added:
- New version of Poisson-based fit statistic for minimization ("Poisson 
maximum-likelihood-ratio statistic", via "`--poisson-mlr`"/"`--mlr`" flag) which is always >= 0, 
and can thus be used with Levenberg-Marquardt minimization, unlike the Cash statistic;
actual fit results should be effectively identical to Cash-statistic fits. 
Thanks to David Streich for pointing out this possibility.

- Full bootstrap-resampling output (i.e., the individual best-fit parameters from
each bootstrap iteration) can now be saved to a text file for later analysis,
via the "`--save-bootstrap`" option. (Thanks to David Streich for suggesting this.)

- Extra minimization algorithms from the NLopt library, specified via
the "`--nlopt`" command-line option (basically, this includes all the "local
derivative-free optimization" algorithms from NLopt; see
[http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms](http://ab-
initio.mit.edu/wiki/index.php/NLopt_Algorithms) for more details). The
Nelder-Mead simplex algorithm is of course still available via its usual
flag ("`--nm`"), and is probably the best of all the NLopt algorithms for
general image-fitting.

- Command-line flag "`--fitstat-only`", which is a synonym for "`--chisquare-only`".

### Changed:
- Minor changes to the wording of help text and error messages.

- Intermediate output from Levenberg-Marquardt fitting now uses the term
"fit statistic" instead of "chi^2" (since L-M can now minimize modified Cash
statistic as well).

- Updates to documentation.



## 1.0.3 -- 2014-07-21
### Fixed:
- Fixes a bug in the newly introduced "`--noisy`" printing mode which
would print scrambled parameter values during the fit if one or more of
the parameters were held fixed. (This bug did *not* affect the actual
fitting or the final output values of the best-fit parameters.) Thanks
to Giulia Savorgnan for spotting the bug so quickly!



## 1.0.2 -- 2014-07-17
### Added:
- A new command-line option "`--loud`" (i.e., the opposite of "`--quiet`"),
which causes the L-M and N-M minimizers to print current model parameter
values during the fitting process (once per iteration for L-M and once
per 100 iterations for N-M simplex).

### Changed:
- Pixels with non-finite values in the image to be fitted are now
automatically masked (instead of causing imfit to complain and quit);
thanks to Giulia Savorgnan for suggesting this.

- Pixels in the input and/or error images which would produce non-finite
weight values (in the case of data-based chi^2 fitting) -- e.g., data
values <= 0 after correcting for any previous sky subtraction, for which
1/sqrt(data) produces non-finite values -- are now ignored *if* they are
masked (instead of causing imfit to complain and quit); thanks to
Francicso Carrera for suggesting this.

- The old random-number-generator code in the DE minimizer has been
replaced with same Mersenne Twister algorithm used elsewhere. This is
significant because the old RNG used the same seed each time, so the
sequence of "random numbers" it generated was always the same... This
means that subsequent fits using DE on the same input will no longer
produce exactly the same output parameters

- Fitting small images (e.g., <~ 200 x 200 pixels) on systems with many
(e.g., > 8) cores is now significantly faster; thanks to Andr&eacute;
Luiz de Amorim for investigating this & figuring out how to make it
happen.

- Error messages are now sent to stderr rather than stdout (probably not
noticeably different unless you redirect the output a lot); some slight
changes in wording of error messages.



## 1.0.1 -- 2014-02-21
### Added:
- Added explicit printing of which statistic is being minimized at the
start of the minimization process.

### Changed:
- Added use of DE for bootstrap resampling if Cash statistic is used and
NLopt library wasn't available.

- Minor improvements to the documentation.

### Fixed:
- Fixed an uninitialized boolean-flag bug which was causing OpenMP failures 
(on at least some 32-bit Linux systems).

- Fixed handling of "`--no-openmp`" compilation option in SConstruct file so
that it actually turns off use of OpenMP (thanks to Sergio Pascual for
identifying this problem).

- Fixed handling of "`--no-nlopt`" compilation option (NO_NLOPT preprocessor
definition) in bootstrap-resampling code (thanks to Guillermo Barro for
identifying this problem).

- Fixed compiled Mac version so that it properly includes static library
code for NLopt (thanks to Giulia Savorgnan for spotting this problem).



## 1.0 -- 2013-11-06
Initial public release.
