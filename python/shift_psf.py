#!/usr/bin/env python3

import copy

psfOrig = [-1,   0,    0.1, 0,   -1,
           0,   0.25, 0.5, 0.25,  0,
           0.1, 0.5,  1.0, 0.5,   0.1,
           0,   0.25, 0.5, 0.25,  0,
           -1,  0,    0.1, 0,     -1]

psfSmall = [0.1,   0.2,    0.1,
            0.2,   0.5,    0.2,
            0.1,   0.2,    0.1,]

padBlank = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


def ShiftAndWrapPSF( psfImage, nRows_psf, nCols_psf, blankImage, nRows_pad, nCols_pad ):
	newImage = copy.copy(blankImage)
	
	centerX_psf = nCols_psf // 2
	centerY_psf = nRows_psf // 2
	centerX_pad = nCols_pad // 2
	centerY_pad = nRows_pad // 2
	
	for i in range(nRows_psf):
		for j in range(nCols_psf):
			psfCol = j
			psfRow = i
			pos_in_psf = i*nCols_psf + j
			destCol = (nCols_pad - centerX_psf + psfCol) % nCols_pad
			destRow = (nRows_pad - centerY_psf + psfRow) % nRows_pad
			#print("(%d,%d) --> (%d,%d)" % (psfCol, psfRow, destCol, destRow))
			newImage[destRow*nCols_pad + destCol] = psfImage[pos_in_psf]
	
	return newImage
	

	
def PrintImage( imageVector, nRows, nCols ):
	for i in range(nRows):
		for j in range(nCols):
			pos = i*nCols + j
			outStr = " %f" % imageVector[pos]
			print(outStr,)
		print()
	print("\n")

print("Original PSF image:")
PrintImage(psfOrig, 5, 5)
print("Original blank padded image:")
PrintImage(padBlank, 10, 10)
print()
testImage = ShiftAndWrapPSF(psfOrig, 5, 5, padBlank, 10, 10)
print("Trial version of shifted/wrapped PSF:")
PrintImage(testImage, 10, 10)
testImage_small = ShiftAndWrapPSF(psfSmall, 3, 3, padBlank, 10, 10)
print("Trial version of shifted/wrapped small PSF:")
PrintImage(testImage_small, 10, 10)

