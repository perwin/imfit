#
# Python code for testing ModelObject handling of parameters and function set

# trial parameter vector:
#   [X0,    Y0,    n,   m_e,   r_e, mu_0,   h,    A0,  sigma,  X0,    Y0,    mu_0, h,   n,   m_e,   r_e
p = [512.0, 513.0, 2.0, 15.60, 4.5, 16.00,  25.0, 14.00, 0.5,  525.0, 526.0, 12.5, 5.5, 1.0, 15.90, 2.1]
# function vector
fv = ["Sersic", "exponential", "Gaussian", "exponential", "Sersic"]
paramSizes = [3, 2, 2, 2, 3]
nFunctions = len(fv)

# vector with starting function numbers (0-based) for function sets
setStart = [0, 3]
nSets = len(setStart)

# convert set-start indicies to setStartFlag: vector specifying whether each function
# is start of new function set or not
setStartFlag = nFunctions * [0]
for i in range(nSets):
	n = setStart[i]
	setStartFlag[n] = 1
print "setStartFlag = ", str(setStartFlag)


def CreateModelImage( params ):

	# test to ensure that setStartFlag[0] == 1, so that we're guaranteed to
	# extract x0,y0 values for the first function
	if setStartFlag[0] != 1:
		print "Oops -- setStartFlag[0] is not = 1 (i.e., first function set",
		print " apparently doesn't start at first function, which is not allowed!)"
		return None
	
	print "Full parameter vector = "
	print p
	
	offset = 0
	
	for n in range(nFunctions):
		# calculate new offset number
		if setStartFlag[n] == 1:
			print "Start of new function set ..."
			# OK, we're at the start of a new function set, so we need to skip over
			# two more parameters [(X0,Y0) for new set]
			x0 = p[offset]
			y0 = p[offset + 1]
			offset += 2
		print "Calling setup for function #%d (%s) with (x,y) = (%f,%f)" % (n, fv[n],
					x0, y0)
		print "\toffset = %d:" % offset
		print "\t" + str(p[offset:])
		offset += paramSizes[n]

	# and then we have the looping through the pixels, calling functions' GetValue()
	# method with pixel x,y


CreateModelImage(p)


