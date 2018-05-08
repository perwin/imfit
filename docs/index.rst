.. Imfit documentation master file, created by
   sphinx-quickstart on Tue May  8 10:54:38 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Imfit
=================================

`Imfit <https://www.mpe.mpg.de/~erwin/code/imfit/>`_ is an open-source C++ program for fitting astronomical images
(primarily of galaxies) in FITS format. It is fast, flexible, and
designed to be easily extended with new functions for model image
components.


**Examples of Use:**

Fitting an image:

    ``$ imfit someimage.fits --config model_description.txt``

Fitting a subsection of the image + convolving the model with a Point-Spread-Function image,
using Differential Evolution as the solver:

    ``$ imfit someimage.fits[500:750,600:800] --config model_description.txt --psf psf.fits --de``

Markov Chain Monte Carlo (MCMC) analysis of the same image + model:

    ``$ imfit-mcmc someimage.fits[500:750,600:800] --config model_description.txt --psf psf.fits``

The main documentation for Imfit is in the PDF file 
`imfit_howto.pdf <https://www.mpe.mpg.de/~erwin/resources/imfit/imfit_howto.pdf>`_.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   imfit_tutorial
   config_file_format
   design-and-architecture
   writing_new_image_functions


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
