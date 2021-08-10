#!/bin/bash

python3 -m pip install --upgrade pip
python3 -m pip install -e  .

# add variables to .bashrc
COGEF_BASEDIR=$(pwd)/src
echo 
echo "installation complete"
echo
echo "add the following lines to your .bashrc file:"
echo "export COGEF_BASEDIR=$COGEF_BASEDIR"
echo 
