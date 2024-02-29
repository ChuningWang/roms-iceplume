#!/bin/bash
#set -u


export CALTYPE=greg


##------------------------------------------------------------------------------
## Some calendar tools...
##------------------------------------------------------------------------------

function valid_date 
{
jd=$( ${SCRIPTDIR}/julday.sh $1 $2 $3 $CALTYPE )
${SCRIPTDIR}/caldat.sh $jd $CALTYPE 
}

function makedate
{ 
jd=$( ${SCRIPTDIR}/julday.sh $1 $2 $3 $CALTYPE )
mdy=$( ${SCRIPTDIR}/caldat.sh $jd $CALTYPE )
m=$( echo $mdy | cut -d " " -f 1 )
d=$( echo $mdy | cut -d " " -f 2 )
y=$( echo $mdy | cut -d " " -f 3 )
echo $( printf "%04d\n" $y)$( printf "%02d\n" $m)$( printf "%02d\n" $d) 
}


function sec2hour
{
secs=$1
h=$(( $secs / 3600 ))
m=$(( ( $secs / 60 ) % 60 ))
s=$(( $secs % 60 ))

echo $( printf "%02d:%02d:%02d\n" $h $m $s )
}

##------------------------------------------------------------------------------
## Date of the beginning of the experiment (in $CALTYPE calendar):
##------------------------------------------------------------------------------
#export DATE_BEGIN_EXP=$( makedate $MONTH_BEGIN_EXP $DAY_BEGIN_EXP $YEAR_BEGIN_EXP )

##------------------------------------------------------------------------------
## Date of the beginning of the experiment (in julian calendar (in days)):
##------------------------------------------------------------------------------
export JDAY_BEGIN_EXP=$( ${SCRIPTDIR}/julday.sh ${MONTH_BEGIN_JOB} ${DAY_BEGIN_JOB} ${YEAR_BEGIN_JOB} $CALTYPE )

##------------------------------------------------------------------------------
# Date of the end of the experiment (in $CALTYPE calendar):
##------------------------------------------------------------------------------
#mdy=$( valid_date $MONTH_END_EXP $(( $DAY_BEGIN_EXP + $EXP_DUR_DAY - 1 )) $YEAR_BEGIN_EXP )
#mdy=$( valid_date $MONTH_END_EXP  $(( $DAY_END_EXP - 1 )) $YEAR_END_EXP )
mdy=$( valid_date  $(( $MONTH_BEGIN_JOB + $NBJOB * $JOB_DUR_MTH )) $(( $DAY_BEGIN_JOB + $NBJOB * $JOB_DUR_DAY - 1 )) $YEAR_BEGIN_JOB )
export MONTH_END_EXP=$( echo $mdy | cut -d " " -f 1 )
export DAY_END_EXP=$(   echo $mdy | cut -d " " -f 2 )
export YEAR_END_EXP=$(  echo $mdy | cut -d " " -f 3 )
#
export DATE_END_EXP=$( makedate $MONTH_END_EXP $DAY_END_EXP $YEAR_END_EXP )
export JDAY_END_EXP=$( ${SCRIPTDIR}/julday.sh ${MONTH_END_EXP} ${DAY_END_EXP} ${YEAR_END_EXP} $CALTYPE )
#
#[ $DATE_END_EXP -lt $DATE_BEGIN_EXP ] && echo "ERROR: DATE_END_EXP ($DATE_END_EXP) must be larger than DATE_BEGIN_EXP ($DATE_BEGIN_EXP)... We stop..." && exit

##------------------------------------------------------------------------------
# Date of the beginning of the job (in $CALTYPE calendar):
##------------------------------------------------------------------------------
export DATE_BEGIN_JOB=$( makedate $MONTH_BEGIN_JOB $DAY_BEGIN_JOB $YEAR_BEGIN_JOB )
#[ $DATE_BEGIN_JOB -lt $DATE_BEGIN_EXP ] && echo "ERROR: DATE_BEGIN_JOB ($DATE_BEGIN_JOB) must be larger than DATE_BEGIN_EXP ($DATE_BEGIN_EXP)... We stop..." && exit
[ $DATE_BEGIN_JOB -gt $DATE_END_EXP ] && echo "ERROR: DATE_BEGIN_JOB ($DATE_BEGIN_JOB) must be smaller than DATE_END_EXP ($DATE_END_EXP)... We stop..." && exit

##------------------------------------------------------------------------------
## julian date of the beginning of the job
##------------------------------------------------------------------------------
export JDAY_BEGIN_JOB=$( ${SCRIPTDIR}/julday.sh ${MONTH_BEGIN_JOB} ${DAY_BEGIN_JOB} ${YEAR_BEGIN_JOB} $CALTYPE )

##------------------------------------------------------------------------------
# estimation of the job duration and the end of the job in agreement with 
# the outputs frequency...
##------------------------------------------------------------------------------
#
## year month day of the end of job
mdy=$( valid_date $(( $MONTH_BEGIN_JOB + $JOB_DUR_MTH )) $(( $DAY_BEGIN_JOB + $JOB_DUR_DAY - 1 )) $YEAR_BEGIN_JOB )
export MONTH_END_JOB=$( echo $mdy | cut -d " " -f 1 )
export DAY_END_JOB=$(   echo $mdy | cut -d " " -f 2 )
export YEAR_END_JOB=$(  echo $mdy | cut -d " " -f 3 )
# julian date of the end of the job
JDAY_END_JOB=$( ${SCRIPTDIR}/julday.sh ${MONTH_END_JOB} ${DAY_END_JOB} ${YEAR_END_JOB} $CALTYPE )
# total job duration
export TOTAL_JOB_DUR=$(( $JDAY_END_JOB - JDAY_BEGIN_JOB + 1 ))

##------------------------------------------------------------------------------
# check if that job duration allow to finish on end_exp
##------------------------------------------------------------------------------
CHECK_Y=${YEAR_BEGIN_JOB}
CHECK_M=${MONTH_BEGIN_JOB}
CHECK_D=${DAY_BEGIN_JOB}
CHECK_DATE=$( makedate $CHECK_M $CHECK_D $CHECK_Y )
while [[ ${CHECK_DATE} -le ${DATE_END_EXP} ]]; do
    mdy=$( valid_date $(( $CHECK_M + $JOB_DUR_MTH )) $(( $CHECK_D + $JOB_DUR_DAY -1 )) $CHECK_Y )
    export CHECK_M=$( echo $mdy | cut -d " " -f 1 )
    export CHECK_D=$( echo $mdy | cut -d " " -f 2 )
    export CHECK_Y=$( echo $mdy | cut -d " " -f 3 )
    CHECK_DATE=$( makedate $CHECK_M $CHECK_D $CHECK_Y )
    if [[ ${CHECK_DATE} -eq ${DATE_END_EXP} ]]; then
        break
    elif [[ ${CHECK_DATE} -gt ${DATE_END_EXP} ]]; then
        echo "Please check your JOB_DURATION so the END DATE is perfectly reached"
        exit
    fi
    mdy=$( valid_date $CHECK_M $(( $CHECK_D + 1 )) $CHECK_Y )
    export CHECK_M=$( echo $mdy | cut -d " " -f 1 )
    export CHECK_D=$( echo $mdy | cut -d " " -f 2 )
    export CHECK_Y=$( echo $mdy | cut -d " " -f 3 )
    CHECK_DATE=$( makedate $CHECK_M $CHECK_D $CHECK_Y )
done

##------------------------------------------------------------------------------
# date of the end of the job
##------------------------------------------------------------------------------
export DATE_END_JOB=$( makedate $MONTH_END_JOB $DAY_END_JOB $YEAR_END_JOB )
if [ $DATE_END_JOB -lt $DATE_BEGIN_JOB ] 
    then
    echo "ERROR: DATE_END_JOB ($DATE_END_JOB) must be larger than DATE_BEGIN_JOB ($DATE_BEGIN_JOB)... We stop..." 
    exit
fi
if [ $DATE_END_JOB -gt $DATE_END_EXP ] 
    then 
    echo "ERROR: DATE_END_JOB ($DATE_END_JOB) must be smaller than DATE_END_EXP ($DATE_END_EXP)... We stop..."
    exit
fi

##------------------------------------------------------------------------------
# define dates for next job...
##------------------------------------------------------------------------------
    JDAY_BEGIN_JOBp1=$(( ${JDAY_END_JOB} + 1 ))
    mdy=$( ${SCRIPTDIR}/caldat.sh ${JDAY_BEGIN_JOBp1} ${CALTYPE} )
    MONTH_BEGIN_JOBp1=$( echo ${mdy} | cut -d " " -f 1 )
    DAY_BEGIN_JOBp1=$(   echo ${mdy} | cut -d " " -f 2 )
    YEAR_BEGIN_JOBp1=$(  echo ${mdy} | cut -d " " -f 3 )


