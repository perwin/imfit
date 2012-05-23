#!/usr/bin/env python
#
# Tests to see whether numpy and pyfits are present (returns 1 if they are,
# 0 if at least one is not accessible)

import sys
try:
	import numpy as np
	import pyfits
	librariesPresent = 1
except ImportError:
	librariesPresent = 0

def main():
	print librariesPresent

if __name__ == '__main__':	
	main()
