#!/bin/bash
#
# --------------------------------------------------
#
# Script the mpmd.conf file to be used with srun
#
# --------------------------------------------------
#set -x
#---------------------------------------------------
mod1=$1
mod2=$2
mod3=$3

NBPROC_1=$4
NBPROC_2=$5
NBPROC_3=$6
#---------------------------------------------------
rm -f run_file

mystartproc=0
myoffset=-1
mod1_Str=''
mod2_Str=''
mod3_Str=''

if [ $mod1 != no ] ; then
    myNBPROC_1=$NBPROC_1
    mystartproc=$(( $mystartproc ))
    myendproc=$(( $myNBPROC_1 - 1 ))
    mod1_Str=$mystartproc"-"$myendproc
    echo "$mod1_Str ./$mod1" >> run_file
    myoffset=0
fi

if [ $mod2 != no ] ; then
    myNBPROC_2=$NBPROC_2
    mystartproc=$(( $myoffset + $myendproc + 1))
    myendproc=$(( $mystartproc +  $myNBPROC_2 -1 ))
    mod2_Str=$mystartproc"-"$myendproc
    echo "$mod2_Str ./$mod2" >> run_file
    myoffset=0
fi

if [ $mod3 != no ] ; then
    myNBPROC_3=$NBPROC_3
    mystartproc=$(( $myoffset + $myendproc +  1 ))
    myendproc=$(( $mystartproc +  $myNBPROC_3 - 1 ))
    mod3_Str=$mystartproc"-"$myendproc
    echo "$mod3_Str ./$mod3" >> run_file
fi
