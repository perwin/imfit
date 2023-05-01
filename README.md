<!-- # Imfit -->

<img src="graphics/logo_for_README.png" width=95%>


Imfit
([http://www.mpe.mpg.de/~erwin/code/imfit/](http://www.mpe.mpg.de/~erwin/code/imfit/)) 
is an open-source C++ program for fitting astronomical images
(primarily of galaxies) in FITS format. It is fast, flexible, and designed to be easily
extended with new functions for model image components.

**Examples of Use:**

Fitting an image:

    $ imfit someimage.fits --config model_description.txt

Fitting a subsection of the image + convolving the model with a Point-Spread-Function image,
using Differential Evolution as the solver:

    $ imfit someimage.fits[500:750,600:800] --config model_description.txt --psf psf.fits --de

Markov Chain Monte Carlo (MCMC) analysis of the same image + model:

    $ imfit-mcmc someimage.fits[500:750,600:800] --config model_description.txt --psf psf.fits

<br>

Pre-compiled binaries for Linux and macOS (a.k.a. Mac OS X), along with
documentation and example files, can be found [at the main web
page](http://www.mpe.mpg.de/~erwin/code/imfit/index.html#downloads), as
can a condensed source-code distribution suitable for local compilations
and for making your own modifications. The most recent full source-code release can
be found [here](https://github.com/perwin/imfit/releases/latest).

(Note on pre-compiled binaries: these are compiled as Intel-native binaries. If you
have an Apple-Silicon-based Mac [a.k.a. M1 or M2], they will run under the Rosetta 2
emulation system.)


<!-- 
ArXiv version:
[![arxiv](http://img.shields.io/badge/arXiv-1408.1097-orange.svg?style=flat)](http://arxiv.org/abs/1408.1097)
 -->

**New:** Preliminary/experimental Python wrapper for Imfit: [PyImfit](https://github.com/perwin/pyimfit)


DOI for current release (v1.9.0), archived at Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.594956.svg)](https://doi.org/10.5281/zenodo.594956)


Astronomy Source Code Library entry:

<a href="http://ascl.net/1408.001"><img src="https://img.shields.io/badge/ascl-1408.001-blue.svg?colorB=262255" alt="ascl:1408.001" /></a>

## Documentation

General online documentation is
[here](https://imfit.readthedocs.io).

The main manual in PDF format:
[imfit_howto.pdf](http://www.mpe.mpg.de/~erwin/resources/imfit/imfit_howto.pdf)

A basic [tutorial](http://www.mpe.mpg.de/~erwin/code/imfit/markdown/index.html) is also available.

Finally, the paper describing Imfit can be found here: 
[Erwin 2015](http://adsabs.harvard.edu/abs/2015ApJ...799..226E).





## About the code at this Github site

[![Build Status](https://travis-ci.org/perwin/imfit.svg?branch=master)](https://travis-ci.org/perwin/imfit)
![Imfit build status](https://github.com/perwin/imfit/workflows/Imfit%20CI/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/imfit/badge/?version=latest)](http://imfit.readthedocs.io/en/latest/?badge=latest)



This repository hosts the complete codebase, including regression
and unit tests, miscellaneous extra tests and notes, and assorted
(nonessential) auxiliary programs. If all you want to do is use Imfit to
model images, it might be simpler to grab the precompiled
versions (or the basic source-code tarball) from [the main Imfit
site](http://www.mpe.mpg.de/~erwin/code/imfit/).

But if you're more comfortable forking things on Github, now you can!
(You can also download the latest official source-code release from this site,
[here](https://github.com/perwin/imfit/releases/latest).)

Feel free to open an issue if you find a bug or want to suggest a feature --
or just email me directly [erwin at sigmaxi.net].


## Requirements

Building Imfit requires the following libraries:

* [FFTW](http://www.fftw.org) (version 3)
* [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html) (version 3)
* [GSL](http://www.gnu.org/software/gsl/) (GNU Scientific Library, version 2.0 or later)
* [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt) -- optional, but strongly recommended

Imfit uses [SCons](http://scons.org) as the build tool,
and should build under any version of GCC from 4.8.1 onward, as well as any 
C++-11-aware version of Clang/LLVM that includes OpenMP support.

See the "Getting and Installing Imfit" chapter of the documentation
-- [imfit_howto.pdf](http://www.mpe.mpg.de/~erwin/resources/imfit/imfit_howto.pdf) -- for more details. 


## License

Imfit is licensed under version 3 of the GNU Public License.
