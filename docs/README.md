# Documentation Directory

This directory holds documentation and documentation-related files for Imfit.

## Imfit Manual

The main Imfit manual is in the LaTeX file imfit\_howto.tex, which can
be used to generate a PDF file via pdflatex. An
up-to-date PDF version generated from the LaTeX file (imfit\_howto.pdf)
is also present.

**Note:** By default, the Imfit manual is built using the LaTeX
MinionPro and MnSymbol packages, which assume you have the Minion Pro
fonts and associated LaTeX-compatible files installed on your system. I have found
[this page](http://kieranhealy.org/blog/archives/2012/11/10/installing-minion-pro/)  
useful for the case of getting Minion Pro fonts to work with LaTeX on Mac OS X.
(Or you can just comment out the lines specifying those packages and fall back on
Computer Modern as your typeface....)

(There are also a few pdflatex-specific microtypography features used by
imfit_howto.tex which probably won't work if for some reason you *aren't* using
pdflatex.)


## Files for Sphinx Documentation on Readthedocs.org

Most of the files in this directory (and in the api_ref/ subdirectory)
are meant to be processed with Sphinx (via the ``make html`` command) to
generate pretty HTML documentation (mainly intended to be done on
readthedocs.org, to produce [the main online docs](https://imfit.readthedocs.io).)

Most of the .rst files are generated from the corresponding Markdown files
via the shell script `convert_md_to_rst.sh`.


## Doxygen Files for API Documentation

Other files in this directory are for generating API documentation with Doxygen.
In order to actually generate the files, the doxygen command should be run from the
*parent* directory thus:

`$ doxygen docs/Doxyfile`


## Other Files

imfit_tutorial.md is a Markdown file used for generating the 
online [Tutorial webpage](https://www.mpe.mpg.de/~erwin/code/imfit/markdown/index.html).
