  SUBROUTINE tool_origindate(netcdfid,varid,date_in_sec) 

  IMPLICIT NONE
#include "netcdf.inc"
  INTEGER       :: netcdfid,varid
  CHARACTER*180 :: units
  CHARACTER*19  :: date
  INTEGER       :: lenstr, luni, indst, ierr
  REAL*8        :: tool_datosec,date_in_sec 

  ierr=nf_get_att_text(netcdfid, varid, 'units', units)
  if (ierr .eq. nf_noerr) then
    luni = lenstr(units)
    if (units(1:6).eq.'second')then
    ! Case units as the format 'seconds since YYYY-MM-DD hh:mm:ss'
      indst=15 ! position to start reading date
    elseif (units(1:3).eq.'day') then
    ! Case units as the format 'days since YYYY-MM-DD hh:mm:ss'
      indst=12 ! position to start reading date
    else
     write(*,'(/1x,2A/6x,A/10x,A)') &
               'TOOL_ORIGINDATE ERROR: ', &
               'unknown units of for time variable' ,&
               'Time variable should follow Netcdf CF format: ',&
               '''seconds(days) since YYYY-MM-DD hh:mm:ss'''
      STOP
    endif
  else
    write(*,'(/1x,A/6x,A/10x,A)') &
            'TOOL_ORIGINDATE ERROR: ',&
            'Time variable should have a units attribute following Netcdf CF format: ',&
            '''seconds(days) since YYYY-MM-DD hh:mm:ss'''
    STOP
  endif

  if (luni<indst) then
    write(*,'(/1x,2A/6x,A/10x,A)') &
             'TOOL_ORIGINDATE ERROR: ',&
             'Could not find any origin date in time var units.',&
             'Time variable should follow Netcdf CF format: ',&
             '''seconds(days) since YYYY-MM-DD hh:mm:ss'''
    STOP
  elseif (luni-indst.eq.3) then
    write(*,'(/1x,4A/1x)') &
        'TOOL_ORIGINDATE: ',&
        'Only origin year is specified, suppose it is ', units(indst:luni),&
        '-01-01 00:00:00'
    date=units(indst:luni)//'-01-01 00:00:00'
  elseif (luni-indst.eq.6) then
    write(*,'(/1x,4A/1x)') &
        'TOOL_ORIGINDATE: ',&
        'Only origin year and month are specified, suppose it is ', units(indst:luni),&
        '/01 00:00:00'
    date=units(indst:luni)//'-01 00:00:00'
  elseif (luni-indst.eq.9) then
    write(*,'(/1x,4A/1x)') &
        'TOOL_ORIGINDATE: ',&
        'Only origin year,month,day are specified, suppose it is ', units(indst:luni),&
        ' 00:00:00'
    date=units(indst:luni)//' 00:00:00'
  elseif (luni-indst.eq.12) then
    write(*,'(/1x,4A/1x)') &
        'TOOL_ORIGINDATE: ', &
        'Only origin year,month,day,hour are specified, suppose it is ', units(indst:luni),':00:00'
    date=units(indst:luni)//':00:00'
  elseif (luni-indst.eq.15) then
     write(*,'(/1x,4A/1x)') &
        'TOOL_ORIGINDATE: ', &
        'Only origin year,month,day,hour,minute are specified, suppose it is ', units(indst:luni),':00'   
    date=units(indst:luni)//':00'
  else
    date=units(indst:luni)  
  endif
 
  date_in_sec=tool_datosec(date)

  END SUBROUTINE TOOL_ORIGINDATE
