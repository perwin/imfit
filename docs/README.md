# Documentation Directory

This directory holds documentation and documentation-related files for Imfit.

## Imfit Manual

The main Imfit manual is in the LaTeX file imfit_howto.tex, which can be used to generate a PostScript or PDF file via latex or pdflatex.
An up-to-date PDF version generated from the LaTeX file (imfit_howto.pdf) is also present.

## Doxygen Files for API Documentation

Most of the files in this directory are for generating API documentation with Doxygen.
In order to actually generate the files, the doxygen command should be run from the
*parent* directory thus:

`$ doxygen docs/Doxyfile`

## Other Files
imfit_tutorial.md is a Markdown file used for generating the online [Tutorial webpage](http://www.mpe.mpg.de/~erwin/code/imfit/markdown/index.html);
currently, the figures for that page are not in the repository.
