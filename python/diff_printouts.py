#!/usr/bin/env python

# Python script for comparing two text files (e.g., save output from running imfit),
# with options for how many lines to skip at beginning or end.
# Prints " OK." if they match, or prints trimmed diff if they don't.
#
# Usage:
# $ diff_printouts.py new_file reference_file

from __future__ import print_function

import sys, os, optparse, difflib

DEFAULT_START = 0


def main( argv=None ):

	usageString = "%prog reference_text_file new_text_file\n"
	parser = optparse.OptionParser(usage=usageString, version="%prog ")
	parser.add_option("--skip-first", type="int", dest="skipFirst", default=DEFAULT_START,
				help="skip the first N lines of lines retrieved from  each file")
	parser.add_option("--skip-last", type="int", dest="skipLast", default=0,
				help="skip the last N lines of lines retrieved from each file")
	parser.add_option("--get-last", type="int", dest="getLast", default=0,
				help="retrieve the last N lines of each file (before applying --skip options")

 	(options, args) = parser.parse_args(argv)
 
	# args[0] = name program was called with
	# args[1] = first actual argument, etc.

	referenceTextFile = args[1]
	newTextFile = args[2]
	refLines = open(referenceTextFile).readlines()
	newLines = open(newTextFile).readlines()

	# Apply our equivalent of "tail -n NN"
	if options.getLast > 0:
		refLines = refLines[-options.getLast:]
		newLines = newLines[-options.getLast:]
	
	startIndex = options.skipFirst
	endIndex = None
	# slight trickiness required: to get through end of list uses [start:], while
	# to get through all but last N value uses [start:-N]
	if options.skipLast > 0:
		endIndex = -options.skipLast

	# create a delta (difflib generator) in context_diff form, with 0 lines of context
	diff = difflib.context_diff(refLines[startIndex:endIndex], newLines[startIndex:endIndex],
				fromfile=referenceTextFile, tofile=newTextFile, n=0)
	difflines = list(diff)
	if len(difflines) == 0:
		print(" OK.")
	else:
		print("   \033[1;31m Failed:\033[0m Diff output:")
		for line in difflines:
			print(line, end="")



if __name__ == '__main__':
	
	main(sys.argv)

