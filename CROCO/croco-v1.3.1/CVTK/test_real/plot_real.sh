#!/bin/bash
#set -x

#- dependancies
REQUIRE="matlab pdfcrop gs"
for i in $REQUIRE
do 
 has_it=$(which $i)
 if [ -z $has_it ]; then
 echo -e "\033[1;31m $i NOT available ... We quit \033[0m" && exit 1
 fi
done

b_n=$(basename ${0})
OPTIND=1

x_n='BASIN CANYON EQUATOR INNERSHELF INTERNAL IGW RIVER SEAMOUNT SHELFRONT SOLITON THACKER OVERFLOW UPWELLING VORTEX JET SHOREFACE SANDBAR RIP SWASH TANK GRAV_ADJ ISOLITON KH_INST TIDAL_FLAT DUNE'

x_d=$(dirname $(dirname $PWD))

#- Choice of the options ---
while getopts :hn:d: V
do
  case $V in
        (h) x_h=${OPTARG};
        echo "Usage      : "${b_n} \
            " [-h] [-n EXAMPLE] [-d ROOT_DIR] ";
        echo " -h               : help";       
        echo " -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all";
        echo " -d ROOTDIR       : Root of the git repository, default : same as CVTK";
        echo "";
        exit 0;;
        (n)  x_n=${OPTARG};;
        (d)  x_d=${OPTARG};;
        (:)  echo ${b_n}" : -"${OPTARG}" option : missing value" 1>&2;
        exit 2;;
        (\?) echo ${b_n}" : -"${OPTARG}" option : not supported" 1>&2;
        exit 2;;
  esac
done
shift $(($OPTIND-1));

LIST_EXAMPLE=$x_n
ROOTDIR=$x_d


NUMBER=${@: -1}
i=0
for EXAMPLE in $LIST_EXAMPLE
  do 
    i=${NUMBER:-$((i=$i+1))}
    echo "-------"
    echo $i
    echo "-------"
    example=$EXAMPLE
    [ "$EXAMPLE" != "IGW" ] && example=$(echo $EXAMPLE |tr '[:upper:]' '[:lower:]')
    myscript="plot_${example}"
    sed  "s/makepdf\(.*\)=\(.*\)0\(.*\)/makepdf=1/g" TEST_CASES/${myscript}.m > tmp.txt && mv tmp.txt TEST_CASES/${myscript}.m
     
    \rm $(echo $EXAMPLE |tr '[:upper:]' '[:lower:]')*.pdf
    \rm $(echo $EXAMPLE |tr '[:lower:]' '[:upper:]')*.pdf
    FILE1=$( ls *.pdf )
    matlab -nodesktop  -nosplash -nodisplay -r "addpath ./TEST_CASES; try,${myscript};end;exit" || exit 3
    FILE2=$( comm -3 <( ls *.pdf ) <( echo "$FILE1" ) )

    [[ -z "${FILE2}"  ]] && exit 4

    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFIXEDMEDIA -sPAPERSIZE=a4 -dPSFitPage -sOutputFile=tmp.pdf	${FILE2} 
    \rm ${FILE2}
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dFIXEDMEDIA -sPAPERSIZE=a4 -dPSFitPage -dSubsetFonts=true -dEmbedAllFonts=true -dPDFSETTINGS=/default  \
-sOutputFile=${EXAMPLE}.pdf                          \
-c "<< \
/EndPage   \
   {    \
      2 eq { pop false }    \
      {    gsave   \
            /Times-Roman 10 selectfont \
         0 -30 moveto ( $( cd $ROOTDIR && git describe --all --long) ) show  \
     /Times-Roman 10 selectfont \
         470 -30 moveto ( $(date "+DATE: %d-%m-%Y TIME: %H:%M:%S") ) show \
         /Times-Bold 18 selectfont  \
         -15 860 moveto ( ${EXAMPLE} ) show  \
         grestore    \
   true   \
} ifelse \
}  bind   \
>> setpagedevice"  -c "<<  /BeginPage  {  0.9 0.9 scale 29.75 42.1 translate  } >> setpagedevice"  -f tmp.pdf 

FILE3=${EXAMPLE}.pdf 
[ $i -eq 1 ] && \rm merged.pdf
[ $i -ne 1 ] && FILE3="merged.pdf $FILE3"
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sPAPERSIZE=a4 -dSubsetFonts=true -dEmbedAllFonts=true -dPDFSETTINGS=/default  \
-sOutputFile=merged_tmp.pdf ${FILE3}  
\mv merged_tmp.pdf merged.pdf
\rm ${EXAMPLE}.pdf tmp.pdf

  done

