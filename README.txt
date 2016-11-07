This archive contains imfit, a program for fitting arbitrary sets of 2D model
functions to astronomical images in FITS format. Also included are makeimage, a
program for generating model FITS images (using, for example, the output of
imfit) and imfit-mcmc, which does Markov chain Monte Carlo analysis.  
See http://www.mpe.mpg.de/~erwin/code/imfit/ (and the documentation
file imfit_howto.pdf) for more information.

If you downloaded a binary-only version (e.g., imfit-x.x-macintel.tar.gz, 
imfit-x.x-linux-64.tar.gz, etc.), copy the executable files imfit, imfit-mcmc, 
and makeimage to someplace convenient on your path. See docs/imfit_howto.pdf 
for detailed notes on how to use them. The original source code can be found at 
the URL given above, and also at https://github.com/perwin/imfit/.

If you downloaded the source-code archive (e.g., imfit-x.x-source.tar.gz), see
docs/imfit_howto.pdf for details on compiling imfit, makeimage, and imfit-mcmc.
(As explained there, you will also need to have, at a minimum, the Python-based 
SCons package and the cfitsio and fftw3 libraries in order to compile the programs.
To compile imfit-mcmc you will also need the GNU Scientific Library.)

Sample configuration files and images can be found in the examples/
subdirectory; see the README_examples.txt file there for details.
