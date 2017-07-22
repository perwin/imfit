# Change Log for Imfit

(Formatting and design based on Olivier Lacan's [Keep a CHANGELOG](http://keepachangelog.com/))

## 1.5.0 -- 2017-mm-dd
### Added:

- Multiple oversampling regions: Imfit can now
optionally convolve multiple subsections of the model image with oversampled PSFs,
instead of just one subsection as in previous versions. (This applies to all
the individual programs: imfit, imfit-mcmc, and makeimage.)

By adding multiple invocations of `--overpsf_region x1:x2,y1:y2`, each
region will be convolved with the same oversampled PSF (specified by
`--overpsf <fits-filename>` and `overpsf-scale <int>`). Alternately, you
can also add multiple invocations of `--overpsf <fits-filename>` and
`overpsf-scale <int>` to specify that each oversampled region should be
convolved with its own, separate PSF image.

- New image function: FerrersBar3D. XXX

- Command-line flag `--no-normalize` will tell imfit, imfit-mcmc, and makeimage
to *skip* normalization of any input PSF images (standard or oversampled). xxx

### Changed:

- Binaries for MacOS X (or "macOS") are now compiled with GCC 7 instead of GCC 5.
(But if you're compiling the source code yourself, any version of GCC from 4.8
onward will do.)

- I have started (very slowly) updating some of the code to make use of
C++-11 features. This should have no visible effects, but might start to be
relevant for people hacking around with the code. It also means that future
versions of Imfit will require compilers that support C++-11 (e.g., GCC
versions 4.8 or newer).

- Memory-use estimation now includes the effects of multiple PSF oversampling
regions.

### Fixed:

- In describing the PSF oversampling options, the documentation
(erroneously) gave an example of `--overpsf_region [x1:x2,y1:y2]`, when
the correct format was `--overpsf_region x1:x2,y1:y2`. The documentation
has been updated to recommend the latter format; however, the code has
also been changed so that *both* formats are now acceptable.



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
GCC installed via a package manager (e.g., Homebrew, Fink, MacPorts); or
2) Use Apple's version of clang without OpenMP support.



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
