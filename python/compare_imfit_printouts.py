#!/usr/bin/env python3
#
# A Python script to compare two text files (tail end of output from imfit) and
# allow for small differences in numerical output, as expected for differential
# evolution fits.

from __future__ import print_function


import sys
import optparse


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

def MakeDict( lines ):
	theDict = {}
	nameList = []
	for line in lines:
		pp = line.split()
		if len(pp) > 0:
			theDict[pp[0]] = line
			nameList.append(pp[0])
	return theDict, nameList
	
def RelativeDiff( val1, val2 ):
	return abs((val1 - val2) / val1)

def CompareResultsEqual( textFile1, textFile2 ):
	
	textLines1 = open(textFile1).readlines()
	textLines2 = open(textFile2).readlines()
	if len(textLines1) != len(textLines1):
		return False

	dict1, names1 = MakeDict(textLines1)
	dict2, names2 = MakeDict(textLines2)
	
	for name in names1:
		if name not in names2:
			print("\tLine beginning with \"{0}\" is missing from file {1}!".format(name, textFile2))
			return False
		line1 = dict1[name]
		line2 = dict2[name]
		if name in ['Reduced', 'AIC', 'X0', 'Y0', 'PA', 'ell', 'n', 'I_e', 'r_e']:
			#proper line with numerical values
			if name in ['X0', 'Y0', 'PA', 'ell', 'n', 'I_e', 'r_e']:
				# parameter lines; we use tolerance = 5e-5 to allow for 2--3e-5 differences in PA
				num1,num2 = float(line1.split()[1]), float(line2.split()[1])
				tolerance = 5.0e-5
				relDiff = RelativeDiff(num1,num2)
				if (relDiff > tolerance):
					print("\tValue of {0} differs by {1}".format(name, tolerance))
					return False
			else:
				tolerance = 1.0e-9
				if name == "Reduced":   # Reduced chi^2 line
					num1,num2 = float(line1.split()[3]), float(line2.split()[3])
				else:   # AIC, BIC line
					n1 = line1.split()[2].rstrip(",")
					n2 = line2.split()[2].rstrip(",")
					num1,num2 = float(n1), float(n2)
				relDiff = RelativeDiff(num1,num2)
				if (relDiff > tolerance):
					print("\tValue of {0} differs by {1}".format(name, tolerance))
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
