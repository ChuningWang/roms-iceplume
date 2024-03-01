#!/bin/bash

cd ${SCRIPTDIR}/../PREPRO/XIOS

if [[ ${USE_XIOS_OCE} == 1 ]]; then
     sed -e "s/OCE_XIOS=.*/OCE_XIOS=\"TRUE\"/g" process_xios_xml.sh > process_xios_xml.tmp
else
     sed -e "s/OCE_XIOS=.*/OCE_XIOS=\"FALSE\"/g" process_xios_xml.sh > process_xios_xml.tmp
fi
mv process_xios_xml.tmp process_xios_xml.sh
chmod 755 process_xios_xml.sh

if [[ ${USE_CPL} -ge 1 ]]; then
    sed -e "s/USE_OASIS=.*/USE_OASIS=\"TRUE\"/g" process_xios_xml.sh > process_xios_xml.tmp
else
    sed -e "s/USE_OASIS=.*/USE_OASIS=\"FALSE\"/g" process_xios_xml.sh > process_xios_xml.tmp
fi
mv process_xios_xml.tmp process_xios_xml.sh
chmod 755 process_xios_xml.sh

if [[ ${USE_XIOS_ATM} == 1 ]]; then
    sed -e "s/ATM_XIOS=.*/ATM_XIOS=\"TRUE\"/g" process_xios_xml.sh > process_xios_xml.tmp
else
    sed -e "s/ATM_XIOS=.*/ATM_XIOS=\"FALSE\"/g" process_xios_xml.sh > process_xios_xml.tmp
fi
mv process_xios_xml.tmp process_xios_xml.sh
chmod 755 process_xios_xml.sh

./process_xios_xml.sh >& log.process_xml

cd -
