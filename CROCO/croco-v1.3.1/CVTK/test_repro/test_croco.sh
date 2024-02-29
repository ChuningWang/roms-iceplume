#!/bin/bash

###########################
source configure_file
###########################

#set -x 
# Get revision of sources
#-------------------------
numrev0=`sed -n '/revision/{n;p;}' gitinfos`
numrev=`echo $numrev0 | tr -d [:blank:]`
#echo
#echo  "Testing CROCO Rev$numrev"

# Define RVTK
#------------
rm -f Recap_*
touch Recap_${TEST_NAME}.git${numrev}
rm -f mylog.txt
touch mylog.txt

SCRIPT_RVTK=rvtk_fast_qsub.bash
echo -e "   - Run Tests" 
./$SCRIPT_RVTK > Recap_${TEST_NAME}.git${numrev}

# Cleaning : remove the binary check_file
# at the end of the test
rm -Rf check_file_*

#
cat mylog.txt
if [ $? -gt 0 ]; then
    echo "EXITING"
    exit
fi  	

