#!/usr/bin/env python3
#
# A Python script to compare two FITS image files to see if the image data matches
# (assuming that some of the headers might not, so we can't simply do a binary
# file comparison).
# Can also be used to compare the sum of two images to a third, to see if they
# match within some tolerance (currently hard-coded as 10^-6).
#
# The script uses sys.exit() to return either 0 for success or 1 for some kind
# of failure; this is for use with shell scripts for regression tests, etc.

from __future__ import print_function


import sys, os, optparse
import numpy as np
try:
	from astropy.io.fits import open as fits_open
except:
	from pyfits import open as fits_open
#import pyfits

# predefine some ANSI color codes
RED  = '\033[31m' # red
NC = '\033[0m' # No Color


TOLERANCE = 1e-6
TOERLANCE_ALLCLOSE = 1.0e-12


def CompareImagesNew( fname1, fname2, tol=TOERLANCE_ALLCLOSE ):
	imdata1 = fits_open(fname1)[0].data
	imdata2 = fits_open(fname2)[0].data
	return np.allclose(imdata1, imdata2, rtol=0, atol=tol)


# we assume that the original/reference image is the *second* one
def CompareImagesEqual( fname1, fname2, minValue=0.0 ):
	imdata1 = fits_open(fname1)[0].data
	imdata2 = fits_open(fname2)[0].data
	validPix = imdata2 > minValue
	return np.array_equal(imdata1[validPix], imdata2[validPix])


def CompareSum( fname1, fname2, referenceSum_fname, minValue=0.0 ):
	"""Sum the first two images and compare the result with the third; if the maximum
	relative deviation is >= 1e-6, return False, else return True.
	"""
	# kludge to avoid division by zero or very small numbers problem when one or
	# more input image has values equal to or very close to zero: add 1 to each
	# input image (and 2 to the reference image)
	imdata1 = fits_open(fname1)[0].data + 1.0
	imdata2 = fits_open(fname2)[0].data + 1.0
	imSum = imdata1 + imdata2
	refSum_imdata = fits_open(referenceSum_fname)[0].data + 2.0
	devianceImdata = np.abs((imSum / refSum_imdata) - 1.0)
		
	if np.max(devianceImdata) >= TOLERANCE:
# 		i,j = np.unravel_index(devianceImdata.argmax(), devianceImdata.shape)
# 		print()
# 		print(np.max(devianceImdata))
# 		print(imSum[i,j], refSum_imdata[i,j], devianceImdata[i,j])
# 		print()
		return False
	else:
		return True




def main(argv=None):

	usageString = "%prog FITS_file_1 FITS_file_2\n"
	usageString += "OR: %prog --compare-sum FITS_file_1 FITS_file_2 reference_sum_FITS_file\n"
	parser = optparse.OptionParser(usage=usageString, version="%prog ")
	parser.add_option("--compare-sum", action="store_true", dest="compareSum", default=False, help="test that sum of first two images matches third image within tolerances")
	parser.add_option("--min-value", type="float", dest="minValue", default=0.0, help="only test pixels with values > min-value in original image")

	(options, args) = parser.parse_args(argv)
 
	# args[0] = name program was called with
	# args[1] = first actual argument, etc.

	fitsFile1 = args[1]
	fitsFile2 = args[2]
	if not os.path.exists(fitsFile1):
		msg = "unable to find FITS image file %s!\n" % fitsFile1
		print(RED + "ERROR: " + msg + NC)
		sys.exit(1)
	if not os.path.exists(fitsFile2):
		msg = "unable to find FITS image file %s!\n" % fitsFile2
		print(RED + "ERROR: " + msg + NC)
		sys.exit(1)

	if options.compareSum is True:
		refSumFile = args[3]
		if not os.path.exists(refSumFile):
			msg = "unable to find FITS image file %s!\n" % refSumFile
			print(RED + "ERROR: " + msg + NC)
			#print("ERROR: unable to find FITS image file %s!\n" % refSumFile)
			sys.exit(1)
		print("\tComparing sum of %s + %s with %s... " % (fitsFile1, fitsFile2, refSumFile), end="")
		result = CompareSum(fitsFile1, fitsFile2, refSumFile, options.minValue)
		if (result is False):
			txt = "\n\t" + RED + ">>> WARNING:" + NC
			print(txt + " image %s + image %s DOES NOT match %s!" % (fitsFile1, fitsFile2, refSumFile))
			print("t            (one or more pixels differ by > %.1e in relative terms\n" % (TOLERANCE))
			sys.exit(1)
		else:
			print(" OK.")
			sys.exit(0)
	else:
		txt = "\tComparing images %s and %s... " % (fitsFile1, fitsFile2)
		print(txt, end="")
#		result = CompareImagesEqual(fitsFile1, fitsFile2, minValue=options.minValue)
		result = CompareImagesNew(fitsFile1, fitsFile2)
		if (result is False):
			txt = "\n\t" + RED + ">>> WARNING:" + NC
			print(txt + " images %s and %s DO NOT match!\n" % (fitsFile1, fitsFile2))
			sys.exit(1)
		else:
			print(" OK.")
			sys.exit(0)


if __name__ == '__main__':
	
	main(sys.argv)
