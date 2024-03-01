#!/bin/bash

XIOS_NAM_DIR=
ROOT_DIR=
RUNDIR=`pwd`
CPP1="cpp -traditional"

#set -x
set -e

DATE=`date "+%Y%m%d-%H:%M:%S"`
OLD="old_${DATE}"
if [[ "${XIOS_NAM_DIR}" == "${RUNDIR}" ]];then
    echo "Source (RUNDIR) and Destination (XIOS_NAM_DIR) are the same, do not create new repository"
elif [[ -d ${XIOS_NAM_DIR} ]]; then
   mv ${XIOS_NAM_DIR} ${XIOS_NAM_DIR}_${OLD}
   mkdir $XIOS_NAM_DIR
elif [[ ! -d ${XIOS_NAM_DIR} ]]; then
   mkdir $XIOS_NAM_DIR
fi

##copy xios_lauch.file and README_XIOS
[[ -f ${ROOT_DIR}/XIOS/xios_launch.file ]] && { \cp -f ${ROOT_DIR}/XIOS/xios_launch.file $XIOS_NAM_DIR ;} || { echo "Missing xios example file (xios_launch.file), continue..." ; } 

[[ -f ${ROOT_DIR}/XIOS/README_XIOS ]] && { \cp -f ${ROOT_DIR}/XIOS/README_XIOS $XIOS_NAM_DIR ;} || { echo "Missing README file (README_XIOS), continue..." ; }

if [[ -z ${OCE_XIOS+x} ]] || [[ ${OCE_XIOS} == TRUE ]]; then
    echo "Processing xml for OCE model"
    echo "-----------------------------" 
    ##copy cppdefs.h, set_global_definitions.h and cppdefs_dev.h

    if [[ "${XIOS_NAM_DIR}" == "${RUNDIR}" && -f ${RUNDIR}/cppdefs.h ]]; then
        echo "Source (RUNDIR) and Destination (XIOS_NAM_DIR) are the same and cppdefs already exists. Continue..."
    elif [[ -f ${RUNDIR}/cppdefs.h ]]; then
        cp -f ${RUNDIR}/cppdefs.h $XIOS_NAM_DIR 
       
    else
        echo "Missing cppdefs.h files in ${RUNDIR} ..."
        exit 
    fi
       
#    [[ -f ${RUNDIR}/cppdefs.h ]] && { \cp -f ${RUNDIR}/cppdefs.h $XIOS_NAM_DIR ;} || { echo "Missing cppdefs.h files in ${RUNDIR} ..." ; exit ; }

    if [[ -f ${RUNDIR}/set_global_definitions.h && "${XIOS_NAM_DIR}" != "${RUNDIR}" ]]; then
        \cp -f ${RUNDIR}/set_global_definitions.h $XIOS_NAM_DIR/
    else
        \cp -f ${ROOT_DIR}/OCEAN/set_global_definitions.h $XIOS_NAM_DIR/
    fi 
    #-   
    if [[ -f ${RUNDIR}/cppdefs_dev.h && "${XIOS_NAM_DIR}" != "${RUNDIR}" ]]; then    
        \cp -f ${RUNDIR}/cppdefs_dev.h  $XIOS_NAM_DIR/
    else
        \cp -f ${ROOT_DIR}/OCEAN/cppdefs_dev.h $XIOS_NAM_DIR/
    fi
    ##

    ## cpp-process the xml files (except iodef.xml) for croco
    cd $XIOS_NAM_DIR/
    $CPP1 -P -traditional -imacros cppdefs.h  ${ROOT_DIR}/XIOS/field_def_croco.xml_full_withcpp field_def_croco.xml
    $CPP1 -P -traditional -imacros cppdefs.h  ${ROOT_DIR}/XIOS/file_def_croco.xml_full_withcpp file_def_croco.xml
    cd -

    ## copy the xml files (except iodef.xml) for croco and wrf that do not nedd cpp-process
    \cp -f  ${ROOT_DIR}/XIOS/context_croco.xml $XIOS_NAM_DIR
    \cp -f  ${ROOT_DIR}/XIOS/domain_def_croco.xml $XIOS_NAM_DIR
    #
else
    echo "Do not process xml for OCE model"
    echo "---------------------------------"
fi

if [[ ! -z ${ATM_XIOS+x} ]] && [[ ${ATM_XIOS} == "TRUE" ]]; then
    echo "Copy xml for ATM model"
    echo "-----------------------------"
    \cp -f  ${ROOT_DIR}/XIOS/context_wrf.xml $XIOS_NAM_DIR
    \cp -f  ${ROOT_DIR}/XIOS/file_def_wrf.xml $XIOS_NAM_DIR
else 
    echo "Do not use xios for ATM model"
    echo "------------------------------"
fi

## copy the various iodef.xml*, 
if [[ -z ${OCE_XIOS+x} ]] && [[ -z ${ATM_XIOS+x} ]]; then
    \cp -f  ${ROOT_DIR}/XIOS/iodef.xml_croco_xios $XIOS_NAM_DIR/iodef.xml
elif [[ ${OCE_XIOS} == "TRUE" ]] && [[ ${ATM_XIOS} == "FALSE" ]] && [[ ${USE_OASIS} == "FALSE" ]]; then
    \cp -f  ${ROOT_DIR}/XIOS/iodef.xml_croco_xios $XIOS_NAM_DIR/iodef.xml
elif [[ ${OCE_XIOS} == "FALSE" ]] && [[ ${ATM_XIOS} == "TRUE" ]] && [[ ${USE_OASIS} == "FALSE" ]]; then
    \cp -f  ${ROOT_DIR}/XIOS/iodef.xml_wrf_xios $XIOS_NAM_DIR/iodef.xml
elif [[ ${OCE_XIOS} == "FALSE" ]] && [[ ${ATM_XIOS} == "TRUE" ]] && [[ ${USE_OASIS} == "TRUE" ]]; then
    \cp -f  ${ROOT_DIR}/XIOS/iodef.xml_croco_noxios_wrf_xios $XIOS_NAM_DIR/iodef.xml
elif [[ ${OCE_XIOS} == "TRUE" ]] && [[ ${ATM_XIOS} == "FALSE" ]] && [[ ${USE_OASIS} == "TRUE" ]]; then
    \cp -f  ${ROOT_DIR}/XIOS/iodef.xml_croco_xios_wrf_noxios $XIOS_NAM_DIR/iodef.xml
elif [[ ${OCE_XIOS} == "TRUE" ]] && [[ ${ATM_XIOS} == "TRUE" ]] && [[ ${USE_OASIS} == "TRUE" ]]; then
    \cp -f  ${ROOT_DIR}/XIOS/iodef.xml_croco_xios_wrf_xios $XIOS_NAM_DIR/iodef.xml
else
    echo "ERROR: iodef.xml for your case not found"
fi
