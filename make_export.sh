#!/bin/bash

# Script to bundle up current repo (renamed to imfit-export) and copy it 
# to Vagrant directory for compiling Linux version

linux_vm_dest64="/Users/erwin/vagrant/ubuntu64-16_dev/transfer"
linux_vm_dest32="/Users/erwin/vagrant/ubuntu32-16_dev/transfer"

cd /Users/erwin/coding
rm imfit-export.tar.gz
rm -rf imfit-export
hg clone imfit imfit-export
tar -cf imfit-export.tar imfit-export && gzip imfit-export.tar
cp imfit-export.tar.gz ${linux_vm_dest64}
cp imfit-export.tar.gz ${linux_vm_dest32}
