#!/usr/bin/env python

# script for generating distribution tarballs

import sys, os, shutil, optparse, tarfile, subprocess, copy

VERSION_STRING = "0.9.7"

os_type = os.uname()[0]   # "Darwin", "Linux", etc.
os_machine_type = os.uname()[4]   # "x86-64", etc.

# basic scons command (specifies OpenMP use and static linking)
scons_string = "scons --openmp --static"

SOURCE_TARFILE = "imfit-%s-source.tar.gz" % VERSION_STRING
if (os_type == "Darwin"):   # OK, we're compiling on Mac OS X
	BINARY_TARFILE = "imfit-%s-macintel.tar.gz" % VERSION_STRING
	# and we can do "fat" compilation (combine 32-bit and 64-bit binaries)
	scons_string += " --fat"
else:
	# assume it's Linux
	if os_machine_type == "x86_64":
		BINARY_TARFILE = "imfit-%s-linux-64.tar.gz" % VERSION_STRING
	else:
		BINARY_TARFILE = "imfit-%s-linux-32.tar.gz" % VERSION_STRING


binary_only_files = """
imfit
makeimage
"""

documentation_files = """
imfit_howto.pdf
imfit_howto.tex
"""
documentationFileDict = {"dir": "docs", "file_list": documentation_files.split()}

misc_required_files = """
COPYING.txt
DISCLAIMER
README.txt
"""

source_header_files = """
definitions
function_object
mp_enorm
param_struct
statistics
utilities_pub
mpfit_cpp
DESolver
add_functions
commandline_parser
config_file_parser
convolver
diff_evoln_fit
image_io
model_object
print_results
"""

source_files_c = """
mp_enorm
statistics
"""

# the following are C++ files
source_files = """
model_object
convolver
commandline_parser 
utilities 
image_io 
mpfit 
diff_evoln_fit
DESolver
config_file_parser 
add_functions 
print_results 
imfit_main
makeimage_main
"""

source_files_funcobj = """
function_object 
func_gaussian 
func_exp 
func_gen-exp  
func_sersic 
func_gen-sersic 
func_core-sersic 
func_broken-exp 
func_broken-exp2d
func_moffat 
func_flatsky 
func_gaussian-ring 
func_gaussian-ring2side
func_edge-on-disk_n4762 
func_edge-on-disk_n4762v2 
func_edge-on-ring 
func_edge-on-ring2side
func_edge-on-disk
func_expdisk3d 
func_gaussianring3d 
integrator
"""
funcObjFileDict = {"dir": "function_objects", "file_list": source_files_funcobj.split()}

example_files = """
config_sersic_ic3478_256.dat
ic3478rss_256.fits
ic3478rss_256_mask.fits
makeimage_config_moffat_psf_51_for_ic3478rss.dat
psf_moffat_51.fits
README_examples.txt
"""
exampleFileDict = {"dir": "examples", "file_list": example_files.split()}

testing_scripts = """
do_imfit_tests
do_makeimage_tests
compare_fits_files.py
py_startup_test.py
"""

test_files = """
config_imfit_expdisk32.dat
imfit_config_ic3478_64x64.dat
imfit_config_ic3478_64x64b.dat
imfit_config_n3073.dat
config_makeimage_sersictest512_bad1.dat
config_makeimage_sersictest512_bad2.dat
config_makeimage_sersictest512_bad3.dat
config_makeimage_sersictest512_bad4.dat
config_makeimage_sersictest512_bad5.dat
config_imfit_sersictest512_badlimits1.dat
config_imfit_sersictest512_badlimits2.dat
uniform_image32.fits
testimage_expdisk32.fits
ic3478rss_64x64.fits
n3073rss_small.fits
n3073rss_small_mask.fits
biggertest_orig.fits
gensersictest_orig.fits
sersic+exp_orig.fits
gensersictest612_conv_cutout512.fits
imfit_textout1
imfit_textout2
imfit_textout3
imfit_textout3a
imfit_textout4
imfit_textout5
imfit_textout6
imfit_textout7
imfit_textout8
imfit_textout9
config_biggertest_4c.dat
config_makeimage_gensersic512.dat
config_makeimage_sersic+exp512.dat
psf_moffat_35.fits
makeimage_textout1
makeimage_textout2
makeimage_textout3
makeimage_textout4
makeimage_textout5
makeimage_textout6
makeimage_textout7
makeimage_textout8
makeimage_textout9
makeimage_textout10
"""
testFileDict = {"dir": "tests", "file_list": test_files.split()}

# file lists:
binary_only_file_list = binary_only_files.split()
misc_required_files_list = misc_required_files.split()
testing_scripts_list = testing_scripts.split()

header_file_list = [fname + ".h" for fname in source_header_files.split()]
#print header_file_list
c_file_list = [fname + ".c" for fname in source_files_c.split()]
#print c_file_list
cplusplus_file_list = [fname + ".cpp" for fname in source_files.split()]
#print cplusplus_file_list
toplevel_source_list = c_file_list + cplusplus_file_list + header_file_list

documentation_file_list = [ documentationFileDict["dir"] + "/" + fname for fname in documentationFileDict["file_list"] ]

example_file_list = [ exampleFileDict["dir"] + "/" + fname for fname in exampleFileDict["file_list"] ]
test_file_list = [ testFileDict["dir"] + "/" + fname for fname in testFileDict["file_list"] ]

funcobj_file_list_cpp = [ funcObjFileDict["dir"] + "/" + fname + ".cpp" for fname in funcObjFileDict["file_list"] ]
funcobj_file_list_h = [ funcObjFileDict["dir"] + "/" + fname + ".h" for fname in funcObjFileDict["file_list"] ]
funcobj_file_list_h.append(funcObjFileDict["dir"] + "/" + "definitions.h")
funcobj_file_list = funcobj_file_list_h + funcobj_file_list_cpp


allFileLists = [binary_only_file_list, misc_required_files_list, toplevel_source_list, documentation_file_list,
				example_file_list, testing_scripts_list, test_file_list, funcobj_file_list]
subdirs_list = ["docs", "examples", "tests", "function_objects"]



def TrimSConstruct( ):
	lines = open("SConstruct").readlines()
	nLines = len(lines)
	for i in range(nLines):
		if lines[i].startswith("# *** Other programs"):
			lastLineNumber = i
	outf = open("SConstruct_export", 'w')
	for line in lines[0:lastLineNumber]:
		outf.write(line)
	outf.close()


def MakeDistributionDir( ):
	distDir = "imfit-%s/" % VERSION_STRING
	# create distribution subdirectories, if needed
	if not os.path.exists(distDir):
		os.mkdir(distDir)
	for subdir in subdirs_list:
		if not os.path.exists(distDir + subdir):
			os.mkdir(distDir + subdir)
	# copy files to distribution subdirectory
	for fileList in allFileLists:
		for fname in fileList:
			shutil.copy(fname, distDir + fname)
	# copy misc. files requring renaming
	# copy trimmed version of SConstruct
	TrimSConstruct()
	shutil.copy("SConstruct_export", distDir + "SConstruct")


def MakeBinaries( ):
	# Generate appropriate binaries
	print("Calling scons to generate imfit binary...")
	subprocess.check_output(scons_string + " imfit", shell=True)
	print("Calling scons to generate makeimage binary...")
	subprocess.check_output(scons_string + " makeimage", shell=True)
	

def MakeBinaryDist( ):
	distDir = "imfit-%s/" % VERSION_STRING
	final_file_list = binary_only_file_list + misc_required_files_list + documentation_file_list + example_file_list
	
	print("Generating tar file %s..." % BINARY_TARFILE)
	tar = tarfile.open(BINARY_TARFILE, 'w|gz') 
	for fname in final_file_list:
		tar.add(distDir + fname)
	tar.close()


def MakeSourceDist( ):
	distDir = "imfit-%s/" % VERSION_STRING
	final_file_list = example_file_list + misc_required_files_list + documentation_file_list
	final_file_list += header_file_list
	final_file_list += c_file_list
	final_file_list += cplusplus_file_list
	final_file_list += funcobj_file_list
	final_file_list += testing_scripts_list
	final_file_list += test_file_list
	final_file_list.append("SConstruct")
	
	tar = tarfile.open(SOURCE_TARFILE, 'w|gz') 
	for fname in final_file_list:
		tar.add(distDir + fname)
	tar.close()



def main(argv):

	usageString = "%prog [options]\n"
	parser = optparse.OptionParser(usage=usageString, version="%prog ")

	parser.add_option("--source-only", action="store_false", dest="binaryDist", default=True,
					  help="generate *only* the source distribution tarball")
	parser.add_option("--binary-only", action="store_false", dest="sourceDist", default=True,
					  help="generate *only* the binary-only distribution tarball")

	(options, args) = parser.parse_args(argv)
	
	print("\nMaking distribution directory and copying files into it...")
	if options.binaryDist is True:
		print("Generating binary-only distribution (%s)..." % BINARY_TARFILE)
		MakeBinaries()
		MakeDistributionDir()
		MakeBinaryDist()
	if options.sourceDist is True:
		print("Generating source distribution (%s)..." % SOURCE_TARFILE)
		MakeDistributionDir()
		MakeSourceDist()
	
	print("Done!\n")



if __name__ == "__main__":
	main(sys.argv)
