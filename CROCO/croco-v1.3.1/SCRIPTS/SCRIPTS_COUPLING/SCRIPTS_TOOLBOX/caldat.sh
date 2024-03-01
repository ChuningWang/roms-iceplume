#!/bin/bash
set -u
## NAME:
##	caldat.sh
##
## PURPOSE:
##	Return the calendar date and time given julian date.
##	This is the inverse of the function julday.sh
##
## CALLING SEQUENCE:
##	caldat.sh Julian CalendarType
##
## INPUTS:
##	Julian contains the Julian Day Number of the
##	specified calendar date.
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
##              extends the gregorian calendar backward before this date
##
##              30d: 30 days par month calendar (360 days/year)
##
##              noleap: no leap year calendar (365 days/year every years)
##
## OPTIONS:
##       -h: to print this header
##
## OUTPUTS:
##	Month:	Number of the desired month (1 = January, ..., 12 = December).
##
##	Day:	Number of day of the month.
##
##	Year:	Number of the desired year.
##
## RESTRICTIONS:
##	Beware of the accuracy of the computation for large number
##
## MODIFICATION HISTORY:
##	Based on "Numerical Recipies in Fortran 90" Chapter B1 page 1013.
##       Masson Sebastien. October 2002.
##       Masson Sebastien. February 2003:check input parameters
##       Masson Sebastien. February 2003:
##             check heading 0 of input parameters
##             -h option
##       Masson Sebastien, Apr. 2005
##       add calendar noleap. allow year 0 for calendar 30d and noleap
##
# take care of the option
while getopts h V
  do
  case $V in
      h) fname=`which \caldat.sh`
	  sed -n -e "/^##/p" $fname
	  exit;;
      \?) echo $USAGE
      exit;;
  esac
done
shift $(($OPTIND-1));
# Are the inputs well defined?
jd=${1:?Undefined Julian day}
caltype=${2:?Undefined calendar type}
# remove heading 0 of $jd?
if [ $(( $jd )) -eq 0 ]
then
    jd=0
else
    while [ `echo $jd | cut -b 1` -eq 0 ]
      do
      jd=`echo $jd | cut -b 2-`
    done
fi
#
case $caltype in
    julian)    
	[ $jd -lt 0 ] && ja=$(( $jd + 36525*( 1 - $jd / 36525 ) )) || ja=$jd
	
	jb=$(( $ja + 1524 ))
	
	jc=$(( 10*($jb - 2439870) -1221 ))
	jc=$(( (6680*36525 + 10 * $jc)/36525 ))
	
	jd=$(( 365 * $jc + (25 * $jc)/100 ))
	
	je=$(( (10000*( $jb - $jd ))/306001 ))
    
	id=$(( $jb - $jd -(306001 * $je)/10000 ))
	
	mm=$(( $je - 1 ))
	[ $mm -gt 12 ] && mm=$(( $mm -12 ))
	iyyy=$(( $jc - 4715 ))
	[ $mm -gt 2 ] && iyyy=$(( $iyyy - 1 ))
	[ $iyyy -le 0 ] && iyyy=$(( $iyyy - 1 ))
	[ $jd -lt 0 ] && iyyy=$(( $iyyy - 100*(1 - $jd / 36525) ))
	;;
    greg)
	igreg=2299161
	if [ $jd -ge $igreg ] 
	    then
	    jalpha=$(( (100*( $jd - 1867216)-25)/3652425 ))
	    ja=$(( $jd + 1 + $jalpha - $jalpha / 4 ))
	    else
	    ja=$jd
	fi
	${SCRIPTDIR}/caldat.sh $ja julian
	exit
    	;;
    progreg)
	jalpha=$(( (100*( $jd - 1867216)-25)/3652425 ))
	ja=$(( $jd + 1 + $jalpha - $jalpha / 4 ))
	${SCRIPTDIR}/caldat.sh $ja julian
	exit
    	;;
    30d)
	iyyy=$(( $jd/360 + 1 ))
	mm=$(( ($jd%360)/30 + 1 ))
	id=$(( ($jd%360)%30 ))
	while [[ $id -lt 1 ]]; do
	    mm=$(( $mm-1 ))
	    id=$(( $id+30 ))
	done
	while [[ $mm -lt 1 ]]; do
	    iyyy=$(( $iyyy-1 ))
	    mm=$(( $mm+12 ))
	done
	;;
    noleap)
	iyyy=$(( $jd / 365 + 1 ))
	ndays=$(( $jd % 365 ))
	if [ $ndays -eq 0 ] 
	then
	    iyyy=$(( $iyyy-1 ))
	    mm=12
	    id=31
	else if [ $ndays -le 31 ]
	then
	    mm=1
	    id=$ndays
	else if [ $ndays -le 59 ]
	then
	    mm=2
	    id=$(( $ndays - 31 ))
	else if [ $ndays -le 90 ]
	then
	    mm=3
	    id=$(( $ndays - 59 ))
	else if [ $ndays -le 120 ]
	then
	    mm=4
	    id=$(( $ndays - 90 ))
	else if [ $ndays -le 151 ]
	then
	    mm=5
	    id=$(( $ndays - 120 ))
	else if [ $ndays -le 181 ]
	then
	    mm=6
	    id=$(( $ndays - 151 ))
	else if [ $ndays -le 212 ]
	then
	    mm=7
	    id=$(( $ndays - 181 ))
	else if [ $ndays -le 243 ]
	then
	    mm=8
	    id=$(( $ndays - 212 ))
	else if [ $ndays -le 273 ]
	then
	    mm=9
	    id=$(( $ndays - 243 ))
	else if [ $ndays -le 304 ]
	then
	    mm=10
	    id=$(( $ndays - 273 ))
	else if [ $ndays -le 334 ]
	then
	    mm=11
	    id=$(( $ndays - 304 ))
	else
	    mm=12
	    id=$(( $ndays - 334 ))
	fi
	fi
	fi
	fi
	fi
	fi
	fi
	fi
	fi
	fi
	fi
	fi
	;;
    *)
	echo "Unknown Calendar type..."
	exit
	;;
    
esac
echo $mm $id $iyyy


