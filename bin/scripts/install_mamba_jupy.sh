#!/bin/bash
# to run one-off setting of jupyter notebook with mambalics

if [ -d ~/.local/share/jupyter/kernels/mambalics ]; then
 echo "the settings are already existing. Do not run this twice please"
 exit
fi

mkdir -p ~/.local/share/jupyter/kernels/mambalics
ln -s $LiCSARpath/misc/jupykernel.json ~/.local/share/jupyter/kernels/mambalics/kernel.json

#echo "bagr"
#echo "buldozer" > ~/.jupy_installed
