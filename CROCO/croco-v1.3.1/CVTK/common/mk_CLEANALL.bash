#!/bin/bash

for i in `ls Configure_Test/*` ; do
    ii=`echo $i | cut -d/ -f2-`
    echo '=> Clean the dir: '$ii
    rm -Rf $ii
done
