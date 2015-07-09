#!/usr/bin/env python

# script for generating distribution tarballs

import sys, os, shutil, optparse, tarfile, subprocess, copy
import markdown

# where to copy binary (or source) tarballs when completed (specialized for
# Linux virtual machines on MacBook Pro and local web-site directory for Mac
# versions)
LINUX_DEST = "/media/sf_vbox_shared/"
MAC_DEST = "/Users/erwin/Documents/Working/web_site_new/resources/imfit/"
MAC_DEST_BIN = "/Users/erwin/Documents/Working/web_site_new/resources/imfit/binaries/"

MAC_CHANGELOG_MD = "CHANGELOG.md"
MAC_CHANGELOG_DEST = "/Users/erwin/Documents/Working/web_site_new/resources/imfit/CHANGELOG.html"

VERSION_STRING = "1.3b"

os_type = os.uname()[0]   # "Darwin", "Linux", etc.
os_machine_type = os.uname()[4]   # "x86-64", etc.

# basic scons command (specifies use of OpenMP and static linking)
scons_string = "scons --static"

SOURCE_TARFILE = "imfit-%s-source.tar.gz" % VERSION_STRING
if (os_type == "Darwin"):   # OK, we're compiling on Mac OS X
	BINARY_TARFILE = "imfit-%s-macintel.tar.gz" % VERSION_STRING
	BINARY_TARFILE_OLDMAC = "imfit-%s-macintel_10.6-10.7.tar.gz" % VERSION_STRING
	# ensure that we use GCC-5.1 for all compilation
	scons_string += " --use-gcc"
	scons_string_oldmac = scons_string + " --old-mac"
	SOURCE_COPY_DEST_DIR = MAC_DEST
	BINARY_COPY_DEST_DIR = MAC_DEST_BIN
else:
	# assume it's Linux
	if os_machine_type == "x86_64":
		BINARY_TARFILE = "imfit-%s-linux-64.tar.gz" % VERSION_STRING
	else:
		BINARY_TARFILE = "imfit-%s-linux-32.tar.gz" % VERSION_STRING
	BINARY_COPY_DEST_DIR = LINUX_DEST
	SOURCE_COPY_DEST_DIR = LINUX_DEST


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
CHANGELOG.md
"""

# header files in top-level directory
source_header_files = """
definitions
estimate_memory 
mersenne_twister
mp_enorm
param_struct
statistics
utilities_pub
bootstrap_errors
add_functions
commandline_parser
config_file_parser
convolver
image_io
model_object
print_results
oversampled_region
downsample
option_struct_imfit
option_struct_makeimage
sample_configs
"""

# C files in top-level directory 
source_files_c = """
mersenne_twister
mp_enorm
statistics
"""

# the following are C++ files
# C++ files in top-level directory
source_files = """
model_object
bootstrap_errors
convolver
commandline_parser 
estimate_memory 
utilities 
image_io 
config_file_parser
add_functions
print_results
oversampled_region
downsample
imfit_main
makeimage_main
"""

source_files_solvers ="""
levmar_fit
mpfit 
diff_evoln_fit
DESolver
nmsimplex_fit
nlopt_fit
dispatch_solver
solver_results
"""
solversFileDict = {"dir": "solvers", "file_list": source_files_solvers.split()}

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
func_king
func_moffat 
func_flatsky 
func_gaussian-ring 
func_gaussian-ring2side
func_edge-on-disk_n4762 
func_edge-on-disk_n4762v2 
func_edge-on-ring 
func_edge-on-ring2side
func_edge-on-disk
func_brokenexpdisk3d 
func_expdisk3d 
func_gaussianring3d 
integrator
"""
funcObjFileDict = {"dir": "function_objects", "file_list": source_files_funcobj.split()}


example_files = """
config_exponential_ic3478_256.dat
config_sersic_ic3478_256.dat
ic3478rss_256.fits
ic3478rss_256_mask.fits
config_makeimage_moffat_psf.dat
psf_moffat_51.fits
README_examples.txt
"""
exampleFileDict = {"dir": "examples", "file_list": example_files.split()}

testing_scripts = """
do_imfit_tests
do_makeimage_tests
"""

python_files = """
py_startup_test.py
compare_fits_files.py
compare_imfit_printouts.py
diff_printouts.py
"""
pythonFileDict = {"dir": "python", "file_list": python_files.split()}


test_files = """
config_imfit_expdisk32.dat
imfit_config_ic3478_64x64.dat
imfit_config_ic3478_64x64b.dat
imfit_config_ic3478_64x64c.dat
imfit_config_n3073.dat
config_imfit_pgc35772.dat
config_imfit_poisson.dat
config_imfit_flatsky.dat
config_3x3_flatsky.dat
config_imfit_gauss-oversample-test2.dat
config_makeimage_sersictest512_bad1.dat
config_makeimage_sersictest512_bad2.dat
config_makeimage_sersictest512_bad3.dat
config_makeimage_sersictest512_bad4.dat
config_makeimage_sersictest512_bad5.dat
config_imfit_sersictest512_badlimits1.dat
config_imfit_sersictest512_badlimits2.dat
config_imfit_sersictest512_badlimits3.dat
config_imfit_badparamline.dat
uniform_image32.fits
testimage_expdisk32.fits
testimage_poisson_lowsn20.fits
ic3478rss_64x64.fits
n3073rss_small.fits
n3073rss_small_cps.fits
n3073rss_small_mask.fits
pgc35772_continuum.fits
pgc35772_mask.fits
totalmask_64x64.fits
biggertest_orig.fits
gensersictest_orig.fits
sersic+exp_orig.fits
gensersictest612_conv_cutout512.fits
flatsky_128x128.fits
testimage_3x3_nan.fits
testimage_3x3_onezero.fits
mask_for_onezero.fits
oversamp_test4.fits
psf_standard.fits
psf_oversamp.fits
imfit_textout1
imfit_textout2
imfit_textout3
imfit_textout3b
imfit_textout3c_tail
imfit_textout3d
imfit_textout3d2
imfit_textout4
imfit_textout4b
imfit_textout4c
imfit_textout4d
imfit_textout4e
imfit_textout4e2
imfit_textout5a_tail
imfit_textout5b_tail
imfit_textout5c_tail
imfit_textout5d_tail
imfit_textout6
imfit_textout7a
imfit_textout7b
imfit_textout7c
imfit_textout_bad1
imfit_textout_bad2
imfit_textout_bad3
imfit_textout_bad4
imfit_textout_bad5
imfit_textout_bad6
imfit_textout_bad7
imfit_textout_bad8
imfit_textout_badnloptname
config_biggertest_4c.dat
config_makeimage_gensersic512.dat
config_makeimage_sersic+exp512.dat
config_makeimage_sersic+exp512_nosize.dat
psf_moffat_35.fits
makeimage_textout0
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
makeimage_textout11
makeimage_textout12
makeimage_textout13
makeimage_textout14
makeimage_textout15
"""
testFileDict = {"dir": "tests", "file_list": test_files.split()}


# Lists of files:
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
python_file_list = [ pythonFileDict["dir"] + "/" + fname for fname in pythonFileDict["file_list"] ]
test_file_list = [ testFileDict["dir"] + "/" + fname for fname in testFileDict["file_list"] ]

solvers_file_list_cpp = [ solversFileDict["dir"] + "/" + fname + ".cpp" for fname in solversFileDict["file_list"] ]
solvers_file_list_h = [ solversFileDict["dir"] + "/" + fname + ".h" for fname in solversFileDict["file_list"] ]
solvers_file_list = solvers_file_list_h + solvers_file_list_cpp

funcobj_file_list_cpp = [ funcObjFileDict["dir"] + "/" + fname + ".cpp" for fname in funcObjFileDict["file_list"] ]
funcobj_file_list_h = [ funcObjFileDict["dir"] + "/" + fname + ".h" for fname in funcObjFileDict["file_list"] ]
funcobj_file_list_h.append(funcObjFileDict["dir"] + "/" + "definitions.h")
funcobj_file_list = funcobj_file_list_h + funcobj_file_list_cpp


allFileLists = [binary_only_file_list, misc_required_files_list, toplevel_source_list, documentation_file_list,
				example_file_list, python_file_list, testing_scripts_list, test_file_list, solvers_file_list,
				funcobj_file_list]
allFileLists_source = [misc_required_files_list, toplevel_source_list, documentation_file_list,
				example_file_list, python_file_list, testing_scripts_list, test_file_list, solvers_file_list,
				funcobj_file_list]
subdirs_list = ["docs", "examples", "python", "tests", "function_objects", "solvers"]



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


def MakeDistributionDir( mode="binary" ):
	distDir = "imfit-%s/" % VERSION_STRING
	# create distribution subdirectories, if needed
	if not os.path.exists(distDir):
		os.mkdir(distDir)
	for subdir in subdirs_list:
		if not os.path.exists(distDir + subdir):
			os.mkdir(distDir + subdir)
	# copy files to distribution subdirectory
	if (mode == "binary"):
		fileLists = allFileLists
	else:
		fileLists = allFileLists_source
	for fileList in fileLists:
		print fileList
		for fname in fileList:
			print fname
			shutil.copy(fname, distDir + fname)
	# copy misc. files requring renaming
	# copy trimmed version of SConstruct
	TrimSConstruct()
	shutil.copy("SConstruct_export", distDir + "SConstruct")


def MakeFatFile( file1, file2, outputName ):
	cmdText = "lipo -create {0} {1} -output {2}".format(file1, file2, outputName)
	subprocess.check_output(cmdText, shell=True)
	
def MakeFatBinaries( ):
	"""We need this because we now default to using GCC 5.1, which cannot make
	"fat" binaries by default.
	We only do this for the "oldmac" case (Mac OS X 10.6 and 10.7), since
	no one is trying to run 32-bit programs on 10.8 or later.
	"""

	print("   Calling SCons to generate 32-bit imfit binary for Mac OS 10.6/10.7...")
	subprocess.check_output(scons_string_oldmac + " --32bit" + " imfit", shell=True)
	shutil.move("imfit", "imfit32")
	print("   Calling SCons to generate 64-bit imfit binary for Mac OS 10.6/10.7...")
	subprocess.check_output(scons_string_oldmac + " imfit", shell=True)
	shutil.move("imfit", "imfit64")
	print("Merging into combined 32-bit/64-bit binary...")
	MakeFatFile("imfit32", "imfit64", "imfit")

	print("   Calling SCons to generate 32-bit makeimage binary for Mac OS 10.6/10.7...")
	subprocess.check_output(scons_string_oldmac + " --32bit" + " makeimage", shell=True)
	shutil.move("makeimage", "makeimage32")
	print("   Calling SCons to generate 64-bit makeimage binary for Mac OS 10.6/10.7...")
	subprocess.check_output(scons_string_oldmac + " makeimage", shell=True)
	shutil.move("makeimage", "makeimage64")
	print("Merging into combined 32-bit/64-bit binary...")
	MakeFatFile("makeimage32", "makeimage64", "makeimage")
	

def MakeBinaries( mode=None ):
	# Generate appropriate binaries
	if (mode is None):
		# Mac OS 10.8 or newer, or Linux
		print("   Calling SCons to generate imfit binary...")
		subprocess.check_output(scons_string + " imfit", shell=True)
		print("   Calling SCons to generate makeimage binary...")
		subprocess.check_output(scons_string + " makeimage", shell=True)
	else:
		# Mac OS 10.6 or 10.7
		MakeFatBinaries()
	

def MakeBinaryDist( mode=None ):
	distDir = "imfit-%s/" % VERSION_STRING
	final_file_list = binary_only_file_list + misc_required_files_list + documentation_file_list + example_file_list
	
	if (mode is None):
		# Mac OS 10.8 or newer, or Linux
		binaryTarfile = BINARY_TARFILE
	else:
		# Mac OS 10.6 or 10.7
		binaryTarfile = BINARY_TARFILE_OLDMAC
	print("Generating tar file %s..." % binaryTarfile)
	tar = tarfile.open(binaryTarfile, 'w|gz') 
	for fname in final_file_list:
		tar.add(distDir + fname)
	tar.close()

	print("Copying gzipped tar file %s to %s..." % (binaryTarfile, BINARY_COPY_DEST_DIR))
	shutil.copy(binaryTarfile, BINARY_COPY_DEST_DIR)


def MakeSourceDist( ):
	distDir = "imfit-%s/" % VERSION_STRING
	final_file_list = example_file_list + misc_required_files_list + documentation_file_list
	final_file_list += header_file_list
	final_file_list += c_file_list
	final_file_list += cplusplus_file_list
	final_file_list += funcobj_file_list
	final_file_list += solvers_file_list
	final_file_list += python_file_list
	final_file_list += testing_scripts_list
	final_file_list += test_file_list
	final_file_list.append("SConstruct")
	
	tar = tarfile.open(SOURCE_TARFILE, 'w|gz') 
	for fname in final_file_list:
		tar.add(distDir + fname)
	tar.close()
	
	print("Copying gzipped tar file %s to %s..." % (SOURCE_TARFILE, SOURCE_COPY_DEST_DIR))
	shutil.copy(SOURCE_TARFILE, SOURCE_COPY_DEST_DIR)
	print("Generating HTML version of CHANGELOG.md and copying to %s..." % (MAC_CHANGELOG_DEST))
	markdown.markdownFromFile(input=MAC_CHANGELOG_MD, output=MAC_CHANGELOG_DEST)



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
		if (os_type == "Darwin"):
			print("\nGenerating binary-only Mac distribution (%s)..." % BINARY_TARFILE)
			MakeBinaries()
			MakeDistributionDir(mode="binary")
			MakeBinaryDist()
			print("Generating binary-only Mac distribution for 10.6/10.7 (%s)..." % BINARY_TARFILE_OLDMAC)
			MakeBinaries(mode="oldmac")
			#MakeDistributionDir()
			MakeBinaryDist(mode="oldmac")
		else:
			print("\nGenerating binary-only Linux distribution (%s)..." % BINARY_TARFILE)
			MakeBinaries()
			MakeDistributionDir(mode="binary")
			MakeBinaryDist()
	if options.sourceDist is True:
		print("\nGenerating source distribution (%s)..." % SOURCE_TARFILE)
		MakeDistributionDir(mode="source")
		MakeSourceDist()
	
	print("Done!\n")



if __name__ == "__main__":
	main(sys.argv)

