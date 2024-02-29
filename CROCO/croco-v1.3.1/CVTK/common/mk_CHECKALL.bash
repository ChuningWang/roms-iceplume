#!/bin/bash

for i in `ls Configure_Test/*` ; do
    ii=`echo $i | cut -d/ -f2-`
    echo ' '
    echo '=> Check the dir: '$ii
    grep 'DONE' $ii/*.log
    grep 'BUGBIN' $ii/*.log
done
