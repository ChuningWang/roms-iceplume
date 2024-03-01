  FUNCTION tool_sectodat(time)

   !&E---------------------------------------------------------------------
   !&E                 ***  FUNCTION tool_sectodat  ***
   !&E
   !&E ** Purpose : convert a time expressed in seconds elapsed since date_ref
   !&E              (set in para.*****) into a date of the tool_julien calendar,
   !&E              returned as a 19-character string ("yyyy/mm/dd hh:mm:ss")
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : main, step, sflx_rad (meteosat case)
   !&E
   !&E ** External calls : 
   !&E
   !&E ** Used ij-arrays : 
   !&E
   !&E ** Modified variables : tool_sectodat
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  2004-08-20
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used
!  USE parameters
   IMPLICIT NONE
   INTEGER,PARAMETER     :: rlg=8
   REAL(kind=rlg),PARAMETER :: tref = 59958230400_rlg

   !! * Declarations function
   CHARACTER(len=19)          :: tool_sectodat

   !! * Arguments
   !REAL(kind=rlg),INTENT(in)  :: time
   real :: time
   !! * Local declarations
   LOGICAL               :: bissext
   CHARACTER(len=19)     :: date
   INTEGER               :: annee,mois,jour,heure,minute,seconde,tot_jours
   INTEGER               :: nb_siecles,nb_annees,nb_4siecles,nb_4annees
   INTEGER,DIMENSION(12) :: jours_avt_mois
   INTEGER,DIMENSION(365):: mois_from_jour
   REAL(kind=rlg),PARAMETER :: secs_in_minute=60.0_rlg,secs_in_heure=3600.0_rlg,secs_in_jour=86400.0_rlg
   REAL(kind=rlg),PARAMETER :: secs_in_annee=secs_in_jour*365.0_rlg
   REAL(kind=rlg),PARAMETER :: secs_in_siecle=secs_in_jour*(76.0_rlg*365.0_rlg+24.0_rlg*366.0_rlg)
   REAL(kind=rlg),PARAMETER :: secs_in_4annees=secs_in_jour*(3.0_rlg*365.0_rlg+366.0_rlg)
   REAL(kind=rlg),PARAMETER :: secs_in_4siecles=4.0_rlg*secs_in_siecle+secs_in_jour
   REAL(kind=rlg)           :: total_secs

   DATA mois_from_jour /31*01,28*02,31*03,30*04,31*05,30*06,31*07,31*08,30*09,31*10,30*11,31*12/
   DATA jours_avt_mois /0,31,59,90,120,151,181,212,243,273,304,334/

   !!----------------------------------------------------------------------
   !! * Executable part

   total_secs = time+tref
   annee      = 0
   bissext = (MOD(3,2) == 0)

   ! **** on determine l annee recherchee  *****

   ! on verifie que l annee recherchee n est pas l annee d origine (0000)

   IF(total_secs  >=  (secs_in_annee+secs_in_jour)) THEN
     total_secs = total_secs - (secs_in_annee+secs_in_jour)
     nb_4siecles= INT(total_secs/secs_in_4siecles)         ! on determine le siecle
     total_secs = MOD(total_secs,secs_in_4siecles)
     nb_siecles = INT(total_secs/secs_in_siecle)

     ! on fait le test sur le 31 decembre des siecles bissextiles (ie :2000,2400...)

     IF(nb_siecles==4 .AND. total_secs >= secs_in_4siecles-secs_in_jour) nb_siecles = 3
     total_secs = total_secs - nb_siecles*secs_in_siecle
     annee      = 400*nb_4siecles + 100*nb_siecles
     nb_4annees = INT(total_secs/secs_in_4annees)  ! on determine l annee exacte
     total_secs = MOD(total_secs,secs_in_4annees)
     nb_annees  = INT(total_secs/secs_in_annee)

     ! on fait le test sur le 31 decembre des annees bissextiles (annees 0,4,8,12...)

     IF(nb_annees==4 .AND. total_secs >= secs_in_4annees-secs_in_jour) nb_annees=3
     total_secs = total_secs - nb_annees*secs_in_annee
     annee      = annee + 4*nb_4annees + nb_annees + 1

   ENDIF

   ! **** annee bissextile ? (seuls les siecles divisibles par 400 le sont)

   bissext =(MOD(annee,400) == 0 .OR. (MOD(annee,4) == 0 .AND. MOD(annee,100) /= 0))

   ! **** on determine le nombre de jours dans l annee

   tot_jours  = INT(total_secs/secs_in_jour)  ! nombre de jours dans l annee
   total_secs = MOD(total_secs,secs_in_jour)

   ! **** on determine le mois et le jour exact

   IF (bissext .AND. tot_jours >= 59) THEN
     mois = mois_from_jour(tot_jours)
   ELSE
     mois = mois_from_jour(tot_jours+1)
   ENDIF

   IF (bissext .AND. mois>2) THEN
     jour = tot_jours - jours_avt_mois(mois)
   ELSE
     jour = tot_jours - jours_avt_mois(mois) + 1
   ENDIF

   heure      = INT(total_secs/secs_in_heure)  ! **** on calcule l heure exacte
   total_secs = MOD(total_secs,secs_in_heure)
   minute     = INT(total_secs/secs_in_minute) ! **** on calcule les minutes exactes
   total_secs = MOD(total_secs,secs_in_minute)
   seconde = total_secs                        ! **** ce qui reste = les secondes.

   ! on peut alors construire la chaine de caracteres resultante

   WRITE (date,800) annee,mois,jour,heure,minute,seconde
   tool_sectodat = date

800 FORMAT(i4.4,'-',i2.2,'-',i2.2,' ',2(i2.2,':'),i2.2)

  END FUNCTION tool_sectodat
