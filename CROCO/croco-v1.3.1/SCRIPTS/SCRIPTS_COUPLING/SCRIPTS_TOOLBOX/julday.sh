#!/bin/bash
set -u
#
## NAME:
##	julday.sh
##
## PURPOSE:
##	Calculate the Julian Day Number for a given month, day, and year.
##	This is the inverse of the function CALDAT.
##
## CALLING SEQUENCE:
##	julday.sh Month Day Year CalendarType
##
## INPUTS:
##	Month:	Number of the desired month (1 = January, ..., 12 = December).
##
##	Day:	Number of day of the month.
##
##	Year:	Number of the desired year.Year parameters must be valid
##              values from the civil calendar.  Years B.C.E. are represented
##              as negative integers.  Years in the common era are represented
##              as positive integers.  In particular, note that there is no
##              year 0 in the civil calendar.  1 B.C.E. (-1) is followed by
##              1 C.E. (1).
##
##	CalendarType: The calendar type to use for the computation. 
##
##              5 cases are coded:
##
##              julian: for julian calendar
##              leap year only when the year is divisible by 4
##
##              greg: for gregorian calendar
##              leap year when the year is divisible by 4 but not by 100 
##                             or it is divisible by 400
##               => 1900 was not a leap year but 2000 was
##              Gregorian Calander was adopted on Oct. 15, 1582
##              skipping from Oct. 4, 1582 to Oct. 15, 1582
##              before Oct. 4, 1582 the gregorian calendar corresponds
##              to the julian calendar
##
##              progreg: for proleptic gregorian calendar
##              Same as the gregorian calendar after Oct. 15, 1582 and 
##              extend the gregorian calendar backward before this date
##
##              30d: 30 days par month calendar
##
##              noleap: no leap year calendar (365 days/year every years)
##
## OPTIONS:
##       -h: to print this header
##
## OUTPUTS:
##	JULDAY returns the Julian Day Number
##
## RESTRICTIONS:
##	Beware of the accuracy of the computation for large number
##
## MODIFICATION HISTORY:
##	Based on "Numerical Recipies in Fortran 90" Chapter B1 page 1011.
##       Masson Sebastien. October 2002.
##       Masson Sebastien. February 2003:check input parameters
##       Masson Sebastien. February 2003:
##             check heading 0 of input parameters
##             -h option
##       Masson Sebastien, Aug. 2003
##       fix bug for negative and large values of month values
##       eg. julday.sh 349 1 1970 greg
##       Masson Sebastien, Apr. 2005
##       add calendar noleap. allow year 0 for calendar 30d and noleap
##
# take care of the option
while getopts h V
  do
  case $V in
      h) fname=`which \julday.sh`
	  sed -n -e "/^##/p" $fname
	  exit;;
      \?) echo $USAGE
      exit;;
  esac
done
shift $(($OPTIND-1));
# Are the inputs well defined?
mm=${1:?Undefined month}
id=${2:?Undefined day}
iyyy=${3:?Undefined year}
caltype=${4:?Undefined calendar type}
# remove heading 0 in input parameters?
function rmhd0 {
    in=$1
    [[ $( echo $in | cut -b 1 ) == "-" ]] && { neg=1 ;} || { neg=0 ;}
    if [ $(( $in )) -eq 0 ]
	then
	in=0
    else
	while [[ $( echo $in | cut -b 1 ) == "-" || $( echo $in | cut -b 1 ) -eq 0 ]]
	  do
	  in=`echo $in | cut -b 2-`
	done
    fi
    if [[ $neg == 1 ]]
        then
        in=$(( $in * -1 ))
    fi
    echo $in
}
mm=`rmhd0 $mm`
id=`rmhd0 $id`
iyyy=`rmhd0 $iyyy`
#
# force mm between 1 and 12
if [ $mm -lt 0 ]
then
    iyyy=$(( $iyyy + $mm / 12 - 1 ))
    mm=$(( 12 + $mm % 12 ))
fi
if [ $mm -gt 12 ]
then
    iyyy=$(( $iyyy + $mm / 12 ))
    mm=$(( $mm % 12 ))
fi
#
jy=$iyyy
#
case $caltype in
    julian)
	[ $jy -eq 0 ] && echo "there is no year zero in Calendar" && exit
	[ $jy -lt 0 ] && jy=$(( $jy + 1 ))
	[ $mm -le 2 ] && injanfeb=1 || injanfeb=0

	jy=$(( $jy - $injanfeb ))
	jm=$(( $mm + 1 + 12 * $injanfeb ))
	
	tmp=$(( (36525 * $jy)/100 ))
	[ $tmp -lt 0 ] && julday=$(( $tmp - 1 )) || julday=$tmp
	tmp=$(( (306001 * $jm)/10000 ))
	[ $tmp -lt 0 ] && julday=$(( $julday + $tmp-1 )) || julday=$(( $julday + $tmp ))
	julday=$(( $julday + $id + 1720995 ))    
	;;
    greg)
	[ $jy -eq 0 ] && echo "there is no year zero in Calendar" && exit
	julday=`${SCRIPTDIR}/julday.sh $1 $2 $3 julian`
	igreg=$(( 15+31*(10+12*1582) ))
	if [ $(( $id + 31*($mm + 12 * $iyyy) )) -ge $igreg ]
	    then
	    [ $jy -lt 0 ] && jy=$(( $jy + 1 ))
	    [ $mm -le 2 ] && jy=$(( $jy - 1 ))
	    
	    ja=$(( $jy / 100 ))
	    [ $ja -lt 0 ] && $ja=$(( $ja - 1 ))
	    tmp=$(( $ja / 4 ))
	    [ $tmp -lt 0 ] && tmp=$(( $tmp - 1 ))
	    julday=$(( $julday + 2 - $ja + $tmp ))
	fi
	
	;;
    progreg)
	[ $jy -eq 0 ] && echo "there is no year zero in Calendar" && exit
	julday=`${SCRIPTDIR}/julday.sh $1 $2 $3 julian`
	[ $jy -lt 0 ] && jy=$(( $jy + 1 ))
	[ $mm -le 2 ] && jy=$(( $jy - 1 ))
	
	ja=$(( $jy / 100 ))
	[ $ja -lt 0 ] && $ja=$(( $ja - 1 ))
	tmp=$(( $ja / 4 ))
	[ $tmp -lt 0 ] && tmp=$(( $tmp - 1 ))
	julday=$(( $julday + 2 - $ja + $tmp ))
	
	;;
    30d)
	julday=$(( 360*( $jy -1 ) + 30*( $mm - 1 ) + $id ))
    
	;;
    noleap)
	case $mm in
	    1) julday=$(( 365*( $jy -1 ) + $id )) ;;
	    2) julday=$(( 365*( $jy -1 ) + 31 + $id )) ;;
	    3) julday=$(( 365*( $jy -1 ) + 59 + $id )) ;;
	    4) julday=$(( 365*( $jy -1 ) + 90 + $id )) ;;
	    5) julday=$(( 365*( $jy -1 ) + 120 + $id )) ;;
	    6) julday=$(( 365*( $jy -1 ) + 151 + $id )) ;;
	    7) julday=$(( 365*( $jy -1 ) + 181 + $id ))	;;
	    8) julday=$(( 365*( $jy -1 ) + 212 + $id ))	;;
	    9) julday=$(( 365*( $jy -1 ) + 243 + $id )) ;;
	    10) julday=$(( 365*( $jy -1 ) + 273 + $id )) ;;
	    11) julday=$(( 365*( $jy -1 ) + 304 + $id )) ;;
	    12) julday=$(( 365*( $jy -1 ) + 334 + $id )) ;;
	    *)
	echo "Impossible case !!!"
	exit
	;;		
	esac
    
	;;
    *)
	echo "Unknown Calendar type..."
	exit
	;;
	
esac

echo $julday
