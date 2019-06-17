#!/bin/bash

# Script to bundle up current repo (renamed to imfit-export) and copy it 
# to Vagrant directory for compiling Linux version

linux_vm_dest="/Users/erwin/vagrant/ubuntu-16_test"

cd /Users/erwin/coding
rm imfit-export.tar.gz
rm -rf imfit-export
hg clone imfit imfit-export
tar -cf imfit-export.tar imfit-export && gzip imfit-export.tar
cp imfit-export.tar.gz ${linux_vm_dest}
