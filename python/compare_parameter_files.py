#!/usr/bin/env python
#
# Python script for comparing two different best-fit parameter files from imit.

from __future__ import print_function

import sys, os, optparse


def GetParameterFileInfo( paramFile ):
	"""Retrieves data from imfit best-fit parameter file, including fitting
	info.
	Returns a tuple of (headerLines, functionNameLines, paramLines)
	"""
	
	lines = open(paramFile).readlines()
	nLines = len(lines)
	headerLines = [line for line in lines if line[0] == "#"]
	dataLines = [line for line in lines if line[0] != "#" and len(line) > 1]
	startParamsIndex = [i for i in range(nLines) if lines[i][0:2] == "X0"][0]
	functionNameIndices = [i for i in range(nLines) if lines[i][0:8] == "FUNCTION"]
	functionNameLines = [line for line in dataLines if line[0:8] == "FUNCTION"]
	paramLines = [line for line in dataLines if line[0:8] != "FUNCTION"]

	return (headerLines, functionNameLines, paramLines)



def main(argv=None):

	usageString = "%prog [options] file1 file2\n"
	parser = optparse.OptionParser(usage=usageString, version="%prog ")


# 	parser.add_option("-t", "--turnon", action="store_true", dest="turnon",
# 					  default=False, help="set this option to True")
# 	parser.add_option("--dx", type="float", dest="delta_x", default=0.0,
# 						help="floating-point value")
# 	parser.add_option("--xmin", type="int", dest="x_min", default=-1,
# 						help="integer value")
	
	(options, args) = parser.parse_args(argv)

	# args[0] = name program was called with
	# args[1] = first actual argument, etc.
	if (len(args) < 3):
		print("ERROR in %s: user must supply *two* parameter files for comparison!" % args[0])
		return None
	inputExists = True
	paramFile1 = args[1]
	paramFile2 = args[2]
	if (not os.path.exists(paramFile1)):
		print("ERROR in %s: Unable to find input parameter file \"%s\"!" % (args[0],paramFile1))
		inputExists = False
	if (not os.path.exists(paramFile2)):
		print("ERROR in %s: Unable to find input parameter file \"%s\"!" % (args[0],paramFile2))
		inputExists = False
	if inputExists is False:
		return None

	GetParameterFileInfo(paramFile1)


if __name__ == '__main__':
	
	main(sys.argv)
