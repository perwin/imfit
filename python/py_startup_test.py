#!/usr/bin/env python
#
# Tests to see whether numpy and pyfits are present (returns 1 if they are,
# 0 if at least one is not accessible)

numpyPresent = False
astropyPresent = False
pyfitsPresent = False

# check for numpy
try:
	import numpy
except ImportError:
	numpyPresent = False
else:
	numpyPresent = True


# check for FITS-reading modules
try:
	import astropy.fits.io
except ImportError:
	astropyPresent = False
else:
	astropyPresent = True

if not astropyPresent:
	try:
		import pyfits
	except ImportError:
		pyfitsPresent = False
	else:
		pyfitsPresent = True


def main():
	# output to be read by shell script calling this program:
	#    1 = necessary libraries are present
	#    0 = one or more necessary libraries are *not* present
	if (numpyPresent and (astropyPresent or pyfitsPresent)):
		print 1
	else:
		print 0


if __name__ == '__main__':	
	main()
