  SUBROUTINE tool_decompdate(date,jj,mm,aaaa,hh,minu,sec)
   !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE tool_decompdate  ***
   !&E
   !&E ** Purpose : decompose a date entered as a 19-character string
   !&E              with format "yyyy/mm/dd hh:mm:ss"
   !&E              into integers corresponding to the day (dd), month (mm),
   !&E              year (yyyy), hour (hh), minute (minu) and second (sec)
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : meteo, tool_datosec, main, step
   !&E
   !&E ** External calls : 
   !&E
   !&E ** Used ij-arrays : 
   !&E
   !&E ** Modified variables : jj, mm, aaaa, hh, minu, sec
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2004-08-16
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

   IMPLICIT NONE

   !! * Arguments
   CHARACTER(len=19),INTENT(in)   :: date
   INTEGER          ,INTENT(out)  :: jj,mm,aaaa,hh,minu,sec

   !! * Local declarations
   LOGICAL                        :: cont
   CHARACTER(len=2)               :: mois
   CHARACTER(len=2),DIMENSION(12) :: tab_mois=(/'01','02','03','04','05','06',&
                                                '07','08','09','10','11','12'/)
   INTEGER                        :: i

   !!----------------------------------------------------------------------
   !! * Executable part

   READ(date,800)aaaa,mois,jj,hh,minu,sec
   cont = .TRUE.
   i = 1

   DO WHILE(cont.AND. i <= 12)
     IF(tab_mois(i) == mois) THEN
       cont = .FALSE.
     ELSE
       i = i +1
     ENDIF
   ENDDO

   IF(i <= 12) THEN
     mm = i
   ELSE
     PRINT*,'erreur tool_decompdate : impossible de determiner le mois'
     STOP
   ENDIF

800 FORMAT(i4,1x,a2,1x,i2,1x,2(i2,1x),i2)

  END SUBROUTINE tool_decompdate
