#!/usr/bin/env python
#
# A Python script to compare two text files (tail end of output from imfit) and
# allow for small differences in numerical output, as expected for differential
# evolution fits.

from __future__ import print_function


import sys, os, optparse


# Sample of diff between two different runs of DE subset of do_imfit_tests:
# <   CHI-SQUARE = 4555.947657    (4089 DOF)
# ---
# >   CHI-SQUARE = 4555.947654    (4089 DOF)
# 4c4
# < AIC = 4569.975054, BIC = 4614.172020
# ---
# > AIC = 4569.975051, BIC = 4614.172017
# 9,13c9,13
# < PA		18.2617
# < ell		0.235988
# < n		 2.4003
# < I_e		20.0091
# < r_e		60.7617
# ---
# > PA		18.2612
# > ell		0.235989
# > n		2.40029
# > I_e		20.0092
# > r_e		60.7616

# Format of one of the output files (e.g., tests/imfit_textout3c_tail):
sample = """
  CHI-SQUARE = 4555.947657    (4089 DOF)

Reduced Chi^2 = 1.114196
AIC = 4569.975054, BIC = 4614.172020

X0              32.9439
Y0              34.0933
FUNCTION Sersic
PA              18.2617
ell             0.235988
n                2.4003
I_e             20.0091
r_e             60.7617

Saving best-fit parameters in file "bestfit_parameters_imfit.dat"
Done!

"""


def ApproxEqual( val1, val2, tolerance ):
	relativeDiff = (val1 - val2) / val1
	if (relativeDiff) > tolerance:
		return False
	else:
		return True

def CompareResultsEqual( textFile1, textFile2 ):
	
	textLines1 = open(textFile1).readlines()
	textLines2 = open(textFile2).readlines()
	if len(textLines1) != len(textLines1):
		return False

	# get numerical results
	# chisquare line
	p1, p2 = textLines1[0].split(), textLines2[0].split()
	chi21, chi22 = float(p1[2]), float(p2[2])
	if not ApproxEqual(chi21, chi22, 1e-9):
		return False
	# reduced chisquare line
	p1, p2 = textLines1[2].split(), textLines2[2].split()
	chi21, chi22 = float(p1[3]), float(p2[3])
	if not ApproxEqual(chi21, chi22, 1e-9):
		return False
	# AIC,BIC line
	p1, p2 = textLines1[3].split(), textLines2[3].split()
	aic1, aic2 = float(p1[2].rstrip(",")), float(p2[2].rstrip(","))
	if not ApproxEqual(aic1, aic2, 1e-9):
		return False
	bic1, bic2 = float(p1[5]), float(p2[5])
	if not ApproxEqual(bic1, bic2, 1e-9):
		return False
	# parameter lines; we use tol = 5e-5 to allow for 2--3e-5 differences in PA
	for i in [5,6,8,9,10,11,12]:
		p1,p2 = textLines1[i].split(), textLines2[i].split()
		if (p1[0] != p2[0]):
			return False
		if not ApproxEqual(float(p1[1]), float(p2[1]), 5e-5):
			return False

	return True


def main( argv=None ):

 	usageString = "%prog new_text_file reference_text_file\n"
 	parser = optparse.OptionParser(usage=usageString, version="%prog ")

 	(options, args) = parser.parse_args(argv)
 
	# args[0] = name program was called with
	# args[1] = first actual argument, etc.

	newTextFile = args[1]
	referenceTextFile = args[2]

	result = CompareResultsEqual(newTextFile, referenceTextFile)
	if (result is False):
		txt = "\n\t>>> WARNING: comparison of %s and %s " % (newTextFile, referenceTextFile)
		txt += "shows output DOES NOT match to within specified tolerance!\n"
		print(txt)
	else:
		print(" OK.")



if __name__ == '__main__':
	
	main(sys.argv)
