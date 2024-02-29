#!/bin/bash
#
# 1- Adapt the Makefile to your own architecture => Makefile.myarch
# 2- ln -sf Makefile.myarch Makefile
#
# 4- source ../myenv_mypath.sh
source ../myenv_mypath.sh
# 3- Clean 
rm -Rf toy_atm toy_oce toy_wav
make clean
# 5- Compile the toy_exe
make
# 6- Link toy_exe to toy_atm, toy_oce and toy_wav
ln -sf toy_model toy_atm
ln -sf toy_model toy_oce
ln -sf toy_model toy_wav
