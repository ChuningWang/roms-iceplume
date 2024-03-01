  SUBROUTINE tool_datetosec(jour,mois,annee,heure,minute,seconde, datetosec)

   !&E---------------------------------------------------------------------
   !&E                 ***  SUBROUTINE tool_datetosec  ***
   !&E
   !&E ** Purpose : return the number of seconds elapsed since the origin date date_ref
   !&E
   !&E ** Description : variant of function tool_datetosec for which the date argument
   !&E    is already provided in terms of integers holding the year, months, day,...
   !&E    etc. As such it doesn't need to call tool_decompdate
   !&E
   !&E ** Called by : init_oa subroutine in the module_oa_interface (BLXD 2020)
   !&E
   !&E ** External calls : None
   !&E
   !&E ** Used ij-arrays : 
   !&E
   !&E ** Modified variables :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2004-08-16
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
!  USE parameters
!  USE comvars2d
   IMPLICIT NONE
   INTEGER,PARAMETER     :: rlg=8
   REAL(kind=rlg),PARAMETER :: tref = 59958230400_rlg

   !! * Declaration function
   REAL(kind=rlg), INTENT(OUT)  :: datetosec

   !! * Arguments
   ! CHARACTER(len=19),INTENT(in) :: date

   !! * Local declarations
   INTEGER, INTENT(IN)        :: annee, mois, jour, heure, minute, seconde
   INTEGER,DIMENSION(12),PARAMETER :: jours_avt_mois=(/0,31,59,90,120,151,181,212,243,273,304,334/)
   REAL(kind=rlg),PARAMETER   :: secs_in_minute=60.0_rlg,secs_in_heure=3600.0_rlg
   REAL(kind=rlg),PARAMETER   :: secs_in_jour=86400.0_rlg,secs_in_annee=secs_in_jour*365.0_rlg
   REAL(kind=rlg),PARAMETER   :: secs_in_siecle=secs_in_jour*36524.0_rlg
   REAL(kind=rlg)             :: total_secs

   !!----------------------------------------------------------------------
   !! * Executable part

   !CALL tool_decompdate(date,jour,mois,annee,heure,minute,seconde)

   ! on ajoute le nombre de secondes pour chaque siecle depuis l origine
   total_secs = secs_in_siecle * INT(annee/100)

   ! on ajoute un jour tous les 400 ans (siecle bissextile)
   ! remarque : 0.9975 = (1 - 1/400) pour l annee 0 qui est bissextile
   total_secs = total_secs + secs_in_jour*INT(DBLE(annee)/400.0_rlg+0.9975_rlg)

   ! on ajoute chaque annee depuis le debut du dernier siecle
   total_secs = total_secs + secs_in_annee*MOD(annee,100)

   ! on ajoute un jour pour chaque annee bissextile depuis le dernier siecle
   total_secs = total_secs + secs_in_jour*INT((MOD(annee,100)-1)/4)

   ! on ajoute chaque mois depuis le debut de l annee
   total_secs = total_secs + jours_avt_mois(mois)*secs_in_jour

   ! on ajoute un jour si on est apres fevrier et que l annee est bissextile

   IF(mois > 2) THEN
     IF(MOD(annee,400) == 0) THEN
       total_secs = total_secs + secs_in_jour
     ELSE
       IF ((MOD(annee,4) == 0) .AND. (MOD(annee,100) /= 0)) &
         total_secs = total_secs + secs_in_jour
     ENDIF
   ENDIF

   ! on ajoute le nombre de jours entres en parametre
   total_secs = total_secs + secs_in_jour*(jour-1)

   ! on ajoute le nombre d heures entrees en parametre
   total_secs = total_secs + secs_in_heure*(heure)

   ! on ajoute le nombre de minutes entrees en parametre
   total_secs = total_secs + secs_in_minute*(minute)

   ! on ajoute le nombre de secondes entrees en parametre
   total_secs = total_secs + seconde
   datetosec = total_secs-tref

  END SUBROUTINE tool_datetosec
