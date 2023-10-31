#!/usr/bin/env python3

# script for generating distribution tarballs

import sys, os, shutil, optparse, tarfile, subprocess, copy
import markdown

sys.path.append(os.getcwd())
import distribution_manifest as dm

# where to copy binary (or source) tarballs when completed (specialized for
# Linux virtual machines on MacBook Pro and local web-site directory for Mac
# versions)
#    LINUX_DEST = local folder within Linux VM
LINUX_DEST = "/vagrant/transfer/"
MAC_DEST = "/Users/erwin/Documents/Working/web_site_new/resources/imfit/"
MAC_DEST_BIN = "/Users/erwin/Documents/Working/web_site_new/resources/imfit/binaries/"

MAC_CHANGELOG_MD = "CHANGELOG.md"
MAC_CHANGELOG_DEST = "/Users/erwin/Documents/Working/web_site_new/resources/imfit/CHANGELOG.html"

VERSION_STRING = "1.9.0"

os_type = os.uname()[0]   # "Darwin", "Linux", etc.
os_machine_type = os.uname()[4]   # "x86-64", etc.

# basic scons command (specifies static linking, etc.)
scons_string = "scons --static"
scons_string_mac = "scons --mac-distribution"

SOURCE_TARFILE = "imfit-%s-source.tar.gz" % VERSION_STRING
if (os_type == "Darwin"):   # OK, we're compiling on Mac OS X
	BINARY_TARFILE = "imfit-%s-macintel.tar.gz" % VERSION_STRING
	BINARY_TARFILE_OLDMAC = "imfit-%s-macintel_10.6-10.7.tar.gz" % VERSION_STRING
	scons_string += " --clang-openmp"
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


# Create dictionaries holding subdirectory names and lists of associated files
documentationFileDict = {"dir": "docs", "file_list": dm.documentation_files.split()}
extrasFileDict = {"dir": "extras", "file_list": dm.extras_files.split()}
coreFileDict = {"dir": "core", "file_list": dm.source_files_core.split()}
solversFileDict = {"dir": "solvers", "file_list": dm.source_files_solvers.split()}
mcmcFileDict = {"dir": "cdream", "file_list": dm.source_files_mcmc.split()}
funcObjFileDict = {"dir": "function_objects", "file_list": dm.source_files_funcobj.split()}
exampleFileDict = {"dir": "examples", "file_list": dm.example_files.split()}
pythonFileDict = {"dir": "python", "file_list": dm.python_files.split()}
testFileDict_imfit = {"dir": "tests/imfit_reference", "file_list": dm.test_files_imfit.split()}
testFileDict_mcmc = {"dir": "tests/imfit-mcmc_reference", "file_list": dm.test_files_mcmc.split()}
testFileDict_makeimage = {"dir": "tests/makeimage_reference", "file_list": dm.test_files_makeimage.split()}
testFileDict = {"dir": "tests", "file_list": dm.test_files.split()}


# Lists of files:
binary_only_file_list = dm.binary_only_files.split()
misc_required_files_list = dm.misc_required_files.split()
testing_scripts_list = dm.testing_scripts.split()
python_files_for_binary_dist_list = dm.python_files_for_binary_dist.split()

documentation_file_list = [ documentationFileDict["dir"] + "/" + fname for fname in documentationFileDict["file_list"] ]
extras_file_list = [ extrasFileDict["dir"] + "/" + fname for fname in extrasFileDict["file_list"] ]

example_file_list = [ exampleFileDict["dir"] + "/" + fname for fname in exampleFileDict["file_list"] ]
python_file_list = [ pythonFileDict["dir"] + "/" + fname for fname in pythonFileDict["file_list"] ]
test_file_imfit_list = [ testFileDict_imfit["dir"] + "/" + fname for fname in testFileDict_imfit["file_list"] ]
test_file_mcmc_list = [ testFileDict_mcmc["dir"] + "/" + fname for fname in testFileDict_mcmc["file_list"] ]
test_file_makeimage_list = [ testFileDict_makeimage["dir"] + "/" + fname for fname in testFileDict_makeimage["file_list"] ]
test_file_list = [ testFileDict["dir"] + "/" + fname for fname in testFileDict["file_list"] ]

solvers_file_list_cpp = [ solversFileDict["dir"] + "/" + fname + ".cpp" for fname in solversFileDict["file_list"] ]
solvers_file_list_h = [ solversFileDict["dir"] + "/" + fname + ".h" for fname in solversFileDict["file_list"] ]
solvers_file_list = solvers_file_list_h + solvers_file_list_cpp

mcmc_file_list_cpp = [ mcmcFileDict["dir"] + "/" + fname + ".cpp" for fname in mcmcFileDict["file_list"] ]
mcmc_file_list_h = [ mcmcFileDict["dir"] + "/" + fname + ".h" for fname in dm.source_header_files_mcmc.split() ]
mcmc_file_list = mcmc_file_list_h + mcmc_file_list_cpp

core_file_list_cpp = [ coreFileDict["dir"] + "/" + fname + ".cpp" for fname in coreFileDict["file_list"] ]
core_file_list_h = [ coreFileDict["dir"] + "/" + fname + ".h" for fname in dm.source_header_files_core.split() ]
core_file_list = core_file_list_h + core_file_list_cpp

funcobj_file_list_cpp = [ funcObjFileDict["dir"] + "/" + fname + ".cpp" for fname in funcObjFileDict["file_list"] ]
funcobj_file_list_h = [ funcObjFileDict["dir"] + "/" + fname + ".h" for fname in funcObjFileDict["file_list"] ]
funcobj_file_list = funcobj_file_list_h + funcobj_file_list_cpp


allFileLists = [binary_only_file_list, misc_required_files_list, documentation_file_list, extras_file_list,
				example_file_list, python_file_list, testing_scripts_list, test_file_imfit_list, 
				test_file_mcmc_list, test_file_makeimage_list, test_file_list, solvers_file_list,
				mcmc_file_list, core_file_list, funcobj_file_list]
allFileLists_source = [misc_required_files_list, documentation_file_list, extras_file_list,
				example_file_list, python_file_list, testing_scripts_list, test_file_imfit_list, 
				test_file_mcmc_list, test_file_makeimage_list, test_file_list, solvers_file_list,
				mcmc_file_list, core_file_list, funcobj_file_list]
subdirs_list = ["docs", "extras", "examples", "python", "tests", "tests/osx", "tests/linux", 
				"tests/imfit_reference", "tests/imfit-mcmc_reference", "tests/makeimage_reference", 
				"tests/mcmc_data", "function_objects", "solvers", "cdream", "cdream/include", 
				"cdream/include/rng", "core"]



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
		print(fileList)
		for fname in fileList:
			print("%s -> %s" % (fname, distDir + fname))
			shutil.copy(fname, distDir + fname)
	# copy misc. files requring renaming
	# copy trimmed version of SConstruct
	TrimSConstruct()
	shutil.copy("SConstruct_export", distDir + "SConstruct")


# def MakeFatFile( file1, file2, outputName ):
# 	cmdText = "lipo -create {0} {1} -output {2}".format(file1, file2, outputName)
# 	subprocess.check_output(cmdText, shell=True)
# 	
# def MakeFatBinaries( ):
# 	"""We need this because we now default to using GCC 5.1, which cannot make
# 	"fat" binaries by default.
# 	We only do this for the "oldmac" case (Mac OS X 10.6 and 10.7), since
# 	no one is trying to run 32-bit programs on 10.8 or later.
# 	"""
# 
# 	print("   Calling SCons to generate 32-bit imfit binary for Mac OS 10.6/10.7...")
# 	subprocess.check_output(scons_string_oldmac + " --32bit" + " imfit", shell=True)
# 	shutil.move("imfit", "imfit32")
# 	print("   Calling SCons to generate 64-bit imfit binary for Mac OS 10.6/10.7...")
# 	subprocess.check_output(scons_string_oldmac + " imfit", shell=True)
# 	shutil.move("imfit", "imfit64")
# 	print("Merging into combined 32-bit/64-bit binary...")
# 	MakeFatFile("imfit32", "imfit64", "imfit")
# 
# 	print("   Calling SCons to generate 32-bit makeimage binary for Mac OS 10.6/10.7...")
# 	subprocess.check_output(scons_string_oldmac + " --32bit" + " makeimage", shell=True)
# 	shutil.move("makeimage", "makeimage32")
# 	print("   Calling SCons to generate 64-bit makeimage binary for Mac OS 10.6/10.7...")
# 	subprocess.check_output(scons_string_oldmac + " makeimage", shell=True)
# 	shutil.move("makeimage", "makeimage64")
# 	print("Merging into combined 32-bit/64-bit binary...")
# 	MakeFatFile("makeimage32", "makeimage64", "makeimage")

def MakeOldMacBinaries( ):
	"""This is for making 64-bit binaries for Mac OS X 10.6 or 10.7.
	"""

	print("   Calling SCons to generate 64-bit imfit binary for Mac OS 10.6/10.7...")
	subprocess.check_output("scons -c", shell=True)
	subprocess.check_output(scons_string_oldmac + " imfit", shell=True)
	print("   Calling SCons to generate 64-bit makeimage binary for Mac OS 10.6/10.7...")
	subprocess.check_output(scons_string_oldmac + " makeimage", shell=True)
	

# def MakeBinaries( mode=None ):
# 	# Generate appropriate binaries
# 	print("MakeBinaries: mode = %s" % mode)
# 	if (mode == "Darwin") or (mode == "Linux"):
# 		# Mac OS 10.8 or newer, or Linux
# 		if (mode is "Darwin"):
# 			# use "--allstatic" to force static linking of GCC libraries (esp. libgomp)
# 			scons_string_final = scons_string_mac
# 		else:
# 			scons_string_final = scons_string
# 		print("   Calling SCons to generate imfit binary...")
# 		subprocess.check_output(scons_string_final + " imfit", shell=True)
# 		print("   Calling SCons to generate makeimage binary...")
# 		subprocess.check_output(scons_string_final + " makeimage", shell=True)
# 		print("   Calling SCons to generate imfit-mcmc binary...")
# 		subprocess.check_output(scons_string_final + " imfit-mcmc", shell=True)
# 	else:   # mode = "oldmac"
# 		# Mac OS 10.6 or 10.7
# 		MakeOldMacBinaries()

def MakeBinaries( target="mac" ):
	if target == "mac":
		scons_string_final = scons_string_mac
	else:   # Linux
		scons_string_final = scons_string
	print("   Calling SCons to generate imfit binary...")
	subprocess.check_output(scons_string_final + " imfit", shell=True)
	print("   Calling SCons to generate makeimage binary...")
	subprocess.check_output(scons_string_final + " makeimage", shell=True)
	print("   Calling SCons to generate imfit-mcmc binary...")
	subprocess.check_output(scons_string_final + " imfit-mcmc", shell=True)
	

def MakeBinaryDist( mode=None ):
	distDir = "imfit-%s/" % VERSION_STRING
	final_file_list = binary_only_file_list + misc_required_files_list + \
						python_files_for_binary_dist_list + documentation_file_list + \
						extras_file_list + example_file_list
	
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
	final_file_list += extras_file_list
	final_file_list += funcobj_file_list
	final_file_list += solvers_file_list
	final_file_list += mcmc_file_list
	final_file_list += core_file_list
	final_file_list += python_file_list
	final_file_list += testing_scripts_list
	final_file_list += test_file_imfit_list
	final_file_list += test_file_mcmc_list
	final_file_list += test_file_makeimage_list
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
			MakeBinaries("mac")
			MakeDistributionDir(mode="binary")
			MakeBinaryDist()
# 			print("Generating binary-only Mac distribution for 10.6/10.7 (%s)..." % BINARY_TARFILE_OLDMAC)
# 			MakeBinaries(mode="oldmac")
# 			#MakeDistributionDir()
# 			MakeBinaryDist(mode="oldmac")
		else:
			print("\nGenerating binary-only Linux distribution (%s)..." % BINARY_TARFILE)
			MakeBinaries("Linux")
			MakeDistributionDir(mode="binary")
			MakeBinaryDist()
	if options.sourceDist is True:
		print("\nGenerating source distribution (%s)..." % SOURCE_TARFILE)
		MakeDistributionDir(mode="source")
		MakeSourceDist()
	
	print("Done!\n")



if __name__ == "__main__":
	main(sys.argv)

