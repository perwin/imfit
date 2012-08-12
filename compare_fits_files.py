#!/usr/bin/env python
#
# A Python script to compare two FITS image files to see if the image data matches
# (assuming that some of the headers might not, so we can't simply do a binary
# file comparison).
# Can also be used to compare the sum of two images to a third, to see if they
# match within some tolerance (currently hard-coded as 10^-6).

from __future__ import print_function


import sys, os, optparse
import numpy as np
import pyfits


TOLERANCE = 1e-6



def CompareImagesEqual( fname1, fname2 ):
	imdata1 = pyfits.open(fname1)[0].data
	imdata2 = pyfits.open(fname2)[0].data
	return np.array_equal(imdata1, imdata2)


def CompareSum( fname1, fname2, referenceSum_fname ):
	"""Sum the first two images and compare the result with the third; if the maximum
	relative deviation is >= 1e-6, return False, else return True.
	"""
	imdata1 = pyfits.open(fname1)[0].data
	imdata2 = pyfits.open(fname2)[0].data
	imSum = imdata1 + imdata2
	refSum_imdata = pyfits.open(referenceSum_fname)[0].data
	devianceImdata = np.abs((imSum / refSum_imdata) - 1.0)
	if np.max(devianceImdata) >= TOLERANCE:
		return False
	else:
		return True




def main(argv=None):

 	usageString = "%prog FITS_file_1 FITS_file_2\n"
 	usageString = "OR: %prog --compare-sum FITS_file_1 FITS_file_2 reference_sum_FITS_file\n"
 	parser = optparse.OptionParser(usage=usageString, version="%prog ")
	parser.add_option("--compare-sum", action="store_true", dest="compareSum",
					  default=False, help="test that sum of first two images matches third image within tolerances")

 	(options, args) = parser.parse_args(argv)
 
	# args[0] = name program was called with
	# args[1] = first actual argument, etc.

	fitsFile1 = args[1]
	fitsFile2 = args[2]
	if not os.path.exists(fitsFile1):
		print("ERROR: unable to find FITS image file %s!\n" % fitsFile1)
		return None
	if not os.path.exists(fitsFile2):
		print("ERROR: unable to find FITS image file %s!\n" % fitsFile2)
		return None

	if options.compareSum is True:
		refSumFile = args[3]
		if not os.path.exists(refSumFile):
			print("ERROR: unable to find FITS image file %s!\n" % refSumFile)
			return None
		print("\tComparing sum of images %s and %s with %s... " % (fitsFile1, fitsFile2, refSumFile), end="")
		result = CompareSum(fitsFile1, fitsFile2, refSumFile)
		if (result is False):
			print("\n\t>>> WARNING: image %s + image %s DO NOT match %s!\n" % (fitsFile1, fitsFile2, refSumFile))
		else:
			print(" OK.")
	else:
		txt = "\tComparing images %s and %s... " % (fitsFile1, fitsFile2)
		print(txt, end="")
		result = CompareImagesEqual(fitsFile1, fitsFile2)
		if (result is False):
			print("\n\t>>> WARNING: images %s and %s DO NOT match!\n" % (fitsFile1, fitsFile2))
		else:
			print(" OK.")


if __name__ == '__main__':
	
	main(sys.argv)
