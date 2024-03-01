PROGRAM mpp_optimiz_nc
   !!---------------------------------------------------------------------
   !!
   !!                       PROGRAM MPP_OPTIMIZ_NC
   !!                     ***********************
   !!
   !!  PURPOSE :
   !!  ---------
   !!              This program is build to optimize the domain beakdown into
   !!              subdomain for mpp computing.
   !!              Once the grid size, and the land/sea mask is known, it looks
   !!              for all the possibilities within a range of setting parameters
   !!              and determine the optimal.
   !!
   !!              Optimization is done with respect to the maximum number of
   !!              sea processors and to the maximum numbers of procs (jprocx)
   !!                     
   !!              Optional optimization can be performed takink into account
   !!              the maximum available processor memory ppmcal. This is
   !!              activated if jpmen =1
   !!
   !! history:
   !! --------
   !!       original  : 95-12 (Imbard M) for OPA8.1, CLIPPER
   !!       f90       : 03-06 (Molines JM), namelist as input
   !!                 : 05-05 (Molines JM), bathy in ncdf
   !!                 : 18-05 (Benshila R), adaptation for CROCO 
   !!----------------------------------------------------------------------
   !! * modules used
    USE netcdf

    IMPLICIT NONE

    INTEGER ::  jprocx=250   !: maximum number of proc. (Read from namelist)
    INTEGER ::  jpmem=0      !: memory constraint (1) or no constraint (0)
       !                     !  (use 1 with caution as the memory size of 
       !                     !   the code lays on CROCO estimates ...)
       !
    INTEGER ::  &
         jpk =     31,    & !: vertical levels (namelist)
         xi_rho  ,    & !: I-size of the model (namelist)
         eta_rho ,    & !: J-size of the model (namelist)
         Npts    =  2,    & !: number of ghost cells
         numnam  =  4       !: logical unit for the namelist
    NAMELIST /namspace/ jpk, Npts
    NAMELIST /namproc/ jprocx, jpmem

    INTEGER ::  jpnix ,jpnjx  
    REAL(kind=8) :: xlen,ylen

    !
    ! Following variables are used only if jpmem=1
    REAL(KIND=4) ::  ppmpt ,   &
         ppmcal = 225000000., &  !: maximum memory of one processor for a given machine (in 8 byte words)
         ppmin  = 0.4,         & !: minimum ratio to fill the memory
         ppmax  = 0.9            !: maximum ration to fill the memory
    ! Aleph
    !     PARAMETER(ppmcal= 16000000.)
    !Brodie
    !     PARAMETER(ppmcal=250000000.)
    ! Uqbar
    !     PARAMETER(ppmcal=3750000000.)
    ! Zahir
    !     PARAMETER(ppmcal=225000000.)

    CHARACTER(LEN=80) :: cbathy, &       !: File name of the netcdf bathymetry (namelist)
        &                clvar           !: Variable name in netcdf for the bathy to be read
    CHARACTER(LEN=80) :: covdta, cdum
    NAMELIST /namfile/ cbathy, covdta
    NAMELIST /namparam/ ppmcal, ppmin, ppmax
    
    INTEGER :: iumout = 16
    INTEGER :: ji,jj,jni,jnj,jni2,jnj2
    INTEGER :: imoy,isurf,ivide
    INTEGER :: in
    INTEGER :: ipi,ipj
    INTEGER :: inf10,inf30,inf50,iptx,isw
    INTEGER :: iii,iij,iiii,iijj,iimoy,iinf10,iinf30,iinf50
    !
    INTEGER,DIMENSION(:,:),ALLOCATABLE     ::  ippdi, ippdj ,iidom, ijdom
    INTEGER ::  LLm, MMm, NP_XI,NP_ETA
    !
    REAL(KIND=4)                           ::  zmin,zmax,zper,zmem
    REAL(KIND=4)                           ::  zzmin,zzmax,zperx
    REAL(KIND=4),DIMENSION(:,:),ALLOCATABLE  ::  zmask ! xi_rho - eta_rho
    REAL(KIND=4),DIMENSION(:,:),ALLOCATABLE  ::  ztemp ! xi_rho - eta_rho


   ! CDF stuff
    INTEGER :: ncid, ivarid, dimid, istatus
    LOGICAL :: llbon=.FALSE.

    INTEGER :: jjc
    INTEGER :: chunk_size_X, margin_X, chunk_size_E, margin_E
    INTEGER :: Istrmpi, Iendmpi, Jstrmpi, Jendmpi, i_X, j_E
    INTEGER, DIMENSION(:), ALLOCATABLE :: nldi,nlei, nldj,nlej,icount
    INTEGER, DIMENSION(:), ALLOCATABLE :: nleiv, nldiv,nlejv,nldjv
    !
    ! 0. Initialisation
    ! -----------------
    OPEN(numnam,FILE='namelist')
    REWIND(numnam)
    READ(numnam,namspace)

    REWIND(numnam)
    READ(numnam,namfile)

    REWIND(numnam)
    READ(numnam,namparam)

    REWIND(numnam)
    READ(numnam,namproc)

    ! estimated  code size expressed in number of 3D arrays : 
    ! A3d + A2d + arrays in common (here for Benguela ref)
    !!!!! CAUTION  here it's a super raw estimate
    ppmpt = 6.+30./jpk + 85./jpk + 3.*4/jpk + 20 +10 +2*2 +2 +4+6
    
    jpnix = jprocx ; jpnjx= jprocx

    ALLOCATE (ippdi(jpnix,jpnjx), ippdj(jpnix,jpnjx) )
    ALLOCATE (iidom(jpnix,jpnjx), ijdom(jpnix,jpnjx) )
    ALLOCATE (nlei(jprocx), nldi(jprocx) )
    ALLOCATE (nlej(jprocx), nldj(jprocx) )
    ! empty processors
    ALLOCATE (nleiv(jprocx), nldiv(jprocx) )
    ALLOCATE (nlejv(jprocx), nldjv(jprocx) )
    ALLOCATE (ICOUNT(jprocx) ) 

    OPEN(iumout,FILE='processor.layout')
    
    WRITE(iumout,*)
    WRITE(iumout,*) ' optimisation de la partition'
    WRITE(iumout,*) ' ----------------------------'
    WRITE(iumout,*)
    !
    ! * Read cdf mask file
    !
    clvar = 'mask_rho'  

    INQUIRE( FILE=cbathy, EXIST=llbon )
    IF( llbon ) THEN
       istatus=NF90_OPEN(cbathy,NF90_NOWRITE,ncid)
       istatus =NF90_OPEN(cbathy,NF90_NOWRITE,ncid)
       istatus = nf90_inq_dimid(ncid,'xi_rho',dimid)
       istatus = nf90_inquire_dimension(ncid,dimid,len=xi_rho)
       istatus = nf90_inq_dimid(ncid,'eta_rho',dimid)
       istatus = nf90_inquire_dimension(ncid,dimid,len=eta_rho) 
      ALLOCATE (zmask(0:xi_rho-1,0:eta_rho-1))
       ALLOCATE (ztemp(xi_rho,eta_rho))
  !     ALLOCATE (zmask(xi_rho,eta_rho))
       istatus = NF90_INQ_VARID(ncid,clvar,ivarid)
       istatus = NF90_GET_VAR(ncid,ivarid,ztemp)
       istatus = NF90_CLOSE(ncid)
    ELSE
        PRINT *,' File missing : ', TRIM(cbathy)
        STOP
    ENDIF

    !
    zmask(0:xi_rho-1,0:eta_rho-1)=ztemp
    DO jj=0,eta_rho-1
       DO ji=0,xi_rho-1
          zmask(ji,jj)=  MIN(REAL(1.,kind=4),MAX(REAL(0.,kind=4),zmask(ji,jj)))  ! Old vector coding rule ...
       END DO
    END DO
 
  !   DO jj=1,eta_rho
  !     DO ji=1,xi_rho
  !       zmask(ji,jj)=  MIN(REAL(1.,kind=4),MAX(REAL(0.,kind=4),zmask(ji,jj)))  ! Old vector coding rule ...
  !     END DO
  !  END DO
   
    PRINT *,'Number of pts     :', eta_Rho*eta_rho
    PRINT *,'Number of sea pts :', INT(SUM(zmask))
    PRINT *

    !
    !  0. Main loop on all possible combination of processors up to jprocx
    ! --------------------------------------------------------------------
    iii=1 ; iij=1
    iiii=xi_rho ; iijj=eta_rho
    iptx=0
    iimoy=0
    zzmin=0. ; zzmax=0.
    iinf10=0 ; iinf30=0 ; iinf50=0
    zperx=1.
    in=0
   
    LLm = xi_rho -2!TRANSFER(OBC_WEST ,zdumm)-TRANSFER(OBC_EAST ,zdumm)  ! EW boundary or not
    Mmm = eta_rho-2!TRANSFER(OBC_NORTH,zdumm)-TRANSFER(OBC_SOUTH,zdumm)  ! NS boundary or not 

    
    DO jni=1,jpnix
       DO jnj=1,jpnjx        !
        
          !  1. Global characteristics of the jni x jnj decomposition
          ! ---------------------------------------------------------
          !

          ! Limitation of the maximum number of PE's
          IF(jni*jnj > jprocx) goto 1000
          !
          NP_XI=jni
          NP_ETA=jnj        
          chunk_size_X=(LLm+NP_XI-1)/NP_XI      
          chunk_size_E=(MMm+NP_ETA-1)/NP_ETA

          ! we requiere number of interior points > 3*Ghostcells : WHY +1 ???
          ! => to avoid too small domains at the boundaries
          IF (chunk_size_X < 2*Npts .OR. chunk_size_E < 2*Npts) go to 1000

          ipi=chunk_size_X+2*Npts  ! Interior + Ghost cells
          ipj=chunk_size_E+2*Npts   

      !    IF(((ipi+2)/2 - (ipi+1)/2 )> 0)  go to 1000
      !    IF(((ipj+2)/2 - (ipj+1)/2) > 0)  go to 1000

          ! Memory optimization ?
          isw=0
          zmem=ppmpt*ipi*ipj*jpk
          IF(zmem > ppmcal) go to 1000
          IF(jpmem == 1) THEN
             IF(zmem.GT.ppmax*ppmcal.OR.zmem.LT.ppmin*ppmcal) isw=1
          ENDIF
          IF(isw.EQ.1) go to 1000
          in=in+1
          !
          WRITE(iumout,*) '--> number of CPUs ',jni*jnj
          WRITE(iumout,*) ' '
          WRITE(iumout,*) " NP_XI=",jni ," NP_ETA=",jnj
          WRITE(iumout,*) " Lm= ",ipi-2*Npts ," Mm= ",ipj-2*Npts
          zper=(jni*jnj*ipi*ipj)/float(xi_rho*eta_rho)
          WRITE(iumout,*) " ratio Lm*Mm/global domain ",zper
          !
          ivide=0
          imoy=0
          zmin=1.e+20
          zmax=-1.e+20
          inf10=0
          inf30=0
          inf50=0
         
          !  2. Loop on the CPUS : Compute mpi stuff for each given decomposition
          ! -----------------------------------------------------------------------
          !
          xlen=xi_rho ; ylen= eta_rho
          DO jj=1,jnj
             DO ji=1,jni
                j_E=jj-1
                i_X=ji-1
                margin_X=(NP_XI*chunk_size_X-Llm)/2
                istrmpi=1+i_X*chunk_size_X-margin_X
                iendmpi=istrmpi+chunk_size_X-1
                istrmpi=MAX(istrmpi,1)
                iendmpi=MIN(iendmpi,LLm)
                ! 
                margin_E=(NP_ETA*chunk_size_E-MMm)/2
                jstrmpi=1+j_E*chunk_size_E-margin_E
                jendmpi=jstrmpi+chunk_size_E-1             
                jstrmpi=MAX(jstrmpi,1)   
                jendmpi=MIN(jendmpi,Mmm) 

                ! security chack, maybe useless by construction
                !if (margin_X >=chunk_size_X) go to 1000               
                !if (margin_E >=chunk_size_E) go to 1000

                 xlen= min(xlen,real(iendmpi-istrmpi+1))
                 ylen= min(ylen,real(jendmpi-jstrmpi+1))
                 if(xlen<Npts ) go to 1000    
                 if(ylen<Npts ) go to 1000    
                

                ! Check wet points over the entire domain to preserve the MPI communication stencil ???????
                isurf=0
                DO jnj2=Max(jstrmpi-Npts,1),Min(jendmpi+Npts,Mmm)
                   DO  jni2=Max(istrmpi-Npts,1),Min(iendmpi+Npts,LLm)
                      IF(zmask(jni2,jnj2).EQ.1.) isurf=isurf+1
                   END DO
                END DO

                IF(isurf.EQ.0) THEN
                   ivide=ivide+1
                 ELSE
                   imoy=imoy+isurf
                  ENDIF
                zper=float(isurf)/float(ipi*ipj)   ! additional points for ghost cells
                IF(zmin.GT.zper.AND.isurf.NE.0) zmin=zper
                IF(zmax.LT.zper.AND.isurf.NE.0) zmax=zper
                IF(zper.LT.0.1.AND.isurf.NE.0) inf10=inf10+1
                IF(zper.LT.0.3.AND.isurf.NE.0) inf30=inf30+1
                IF(zper.LT.0.5.AND.isurf.NE.0) inf50=inf50+1
                !
                !
                ! 3. End of the loop on the CPUS, print
                ! ------------------------------------------------
                !
             END DO
          END DO
          WRITE(iumout,*) ' number of CPUs       ',jni*jnj
          WRITE(iumout,*) ' number of sea CPUs   ',jni*jnj-ivide
          WRITE(iumout,*) ' number of land CPUs  ',ivide
          WRITE(iumout,*) ' average overhead     ',float(imoy)/float(jni*jnj-ivide)/float(ipi*ipj)
          WRITE(iumout,*) ' minimum overhead     ',zmin
          WRITE(iumout,*) ' maximum overhead     ',zmax
          WRITE(iumout,*) ' nb of p overhead < 10 % ',inf10
          WRITE(iumout,*) ' nb of p         10 < nb < 30 % ',inf30-inf10
          WRITE(iumout,*) ' nb of p         30 < nb < 50 % ',inf50-inf10 -inf30
          WRITE(iumout,*) ' number of integration points   ', (jni*jnj-ivide)*ipi*ipj
          WRITE(iumout,*) ' nbr of additional pts          ', (jni*jnj-ivide)*ipi*ipj-xi_rho*eta_rho
          zper=float((jni*jnj-ivide))*float(ipi*ipj)/float(xi_rho*eta_rho)
          WRITE(iumout,*) ' % sup                          ',zper
          WRITE(iumout,*)

          ! 
          ! 4. Optimum search
          ! -------------------------
          !
          IF(ivide.GT.iptx) THEN
             iii=jni
             iij=jnj
             iiii=ipi
             iijj=ipj
             iptx=ivide
             iimoy=imoy
             zzmin=zmin
             zzmax=zmax
             iinf10=inf10
             iinf30=inf30
             iinf50=inf50
             zperx=zper
          ELSE IF(ivide.EQ.iptx.AND.zperx.LT.zper) THEN
             iii=jni
             iij=jnj
             iiii=ipi
             iijj=ipj
             iimoy=imoy
             zzmin=zmin
             zzmax=zmax
             iinf10=inf10
             iinf30=inf30
             iinf50=inf50
             zperx=zper
          ENDIF
          !
          ! 5. End of loop on all possible decomposition
          ! --------------------------------------------
          !
        1000 continue
       END DO
    END DO

    !
    ! 6. loop on optimal cpus (iii x jjj) for plotting purposes
    ! ---------------------------------------------------------
    !
    jjc=0
    ivide=0
    imoy=0
    DO jj=1,iij
      DO ji=1,iii
          j_E=jj-1
          i_X=ji-1
          NP_XI=iii
          NP_ETA=iij
          chunk_size_X=(LLm+NP_XI-1)/NP_XI      
          margin_X=(NP_XI*chunk_size_X-Llm)/2
          istrmpi=1+i_X*chunk_size_X-margin_X
          iendmpi=istrmpi+chunk_size_X-1
          istrmpi=MAX(istrmpi,1)
          iendmpi=MIN(iendmpi,LLm)
          ! 
          chunk_size_E=(MMm+NP_ETA-1)/NP_ETA
          margin_E=(NP_ETA*chunk_size_E-MMm)/2
          jstrmpi=1+j_E*chunk_size_E-margin_E
          jendmpi=jstrmpi+chunk_size_E-1             
          jstrmpi=MAX(jstrmpi,1)   
          jendmpi=MIN(jendmpi,Mmm) 
          
          ! Check wet points over the entire domain to preserve the MPI communication stencil
          isurf=0
          DO jnj2=Max(jstrmpi-Npts,1),Min(jendmpi+Npts,Mmm)
             DO  jni2=Max(istrmpi-Npts,1),Min(iendmpi+Npts,LLm)
                 IF(zmask(jni2,jnj2).EQ.1.) isurf=isurf+1
             END DO
          END DO

          IF (isurf.EQ.0) THEN
             ivide=ivide+1
             nldiv(ivide)=istrmpi
             nleiv(ivide)=iendmpi
             nldjv(ivide)=jstrmpi
             nlejv(ivide)=jendmpi
          ELSE
             imoy=imoy+isurf
             jjc=jjc+1
             icount(jjc)=isurf 
             nldi(jjc)=istrmpi
             nlei(jjc)=iendmpi
             nldj(jjc)=jstrmpi
             nlej(jjc)=jendmpi
          ENDIF
          if ( iendmpi ==-1) STOP !return
          !
          !
          ! End of the loop on the optimal CPUS
          ! -----------------------------------
          !
      END DO
    END DO

    !
    ! 7. Print the result
    ! -------------------
    !
    IF(in.EQ.0) THEN
       WRITE(iumout,*) ' the choice could not be made '
       WRITE(iumout,*)
       WRITE(iumout,*) ' the max number of CPUs is too small'
       STOP 
    ENDIF
    WRITE(iumout,*) ' optimum choice'
    WRITE(iumout,*) ' =============='
    WRITE(iumout,*) 
    WRITE(iumout,*) '--> Number of CPUs : NNODES = ',iii*iij-iptx
    WRITE(iumout,*) ' '
    WRITE(iumout,*) " NP_XI =",iii ," NP_ETA =",iij
    WRITE(iumout,*) " Lm =",iiii-2*Npts ," Mm =",iijj-2*Npts
    WRITE(iumout,*) 
    WRITE(iumout,*) ' number of sea CPUs      ',iii*iij-iptx
    WRITE(iumout,*) ' number of land CPUs     ',iptx
    WRITE(iumout,*) ' average overhead        ',float(iimoy)/float(iii*iij-iptx)/float(iiii*iijj)
    WRITE(iumout,*) ' minimum overhead        ',zzmin
    WRITE(iumout,*) ' maximum overhead        ',zzmax
    WRITE(iumout,*) ' nb of overhead p. < 10 %         ', iinf10
    WRITE(iumout,*) ' nb of overhead p. 10 < nb < 30 % ', iinf30-iinf10
    WRITE(iumout,*) ' nb de overhead p  30 < nb < 50 % ', iinf50-iinf10 -iinf30
    WRITE(iumout,*) ' number of integration points     ', (iii*iij-iptx)*iiii*iijj
    WRITE(iumout,*) ' number of additionnal pts        ', (iii*iij-iptx)*iiii*iijj-xi_rho*eta_rho
    WRITE(iumout,*) ' % sup                            ', zperx
    WRITE(iumout,*)
    !

    WRITE(*,*) ' optimum choice'
    WRITE(*,*) ' =============='
    WRITE(*,*) 
    WRITE(*,*) '--> Number of CPUs : NNODES = ',iii*iij-iptx
    WRITE(*,*) ' '
    WRITE(*,*) " NP_XI =",iii ," NP_ETA =",iij
    WRITE(*,*) " Lm =",iiii-2*Npts ," Mm =",iijj-2*npts
    WRITE(*,*) 
    WRITE(*,*) ' number of sea CPUs      ',iii*iij-iptx
    WRITE(*,*) ' number of land CPUs     ',iptx
    WRITE(*,*) ' average overhead        ',float(iimoy)/float(iii*iij-iptx)/float(iiii*iijj)
    WRITE(*,*) ' minimum overhead        ',zzmin
    WRITE(*,*) ' maximum overhead        ',zzmax
    WRITE(*,*) ' nb of overhead p. < 10 %         ', iinf10
    WRITE(*,*) ' nb of overhead p. 10 < nb < 30 % ', iinf30-iinf10
    WRITE(*,*) ' nb de overhead p  30 < nb < 50 % ', iinf50-iinf10 -iinf30
    WRITE(*,*) ' number of integration points     ', (iii*iij-iptx)*iiii*iijj
    WRITE(*,*) ' number of additionnal pts        ', (iii*iij-iptx)*iiii*iijj-xi_rho*eta_rho
    WRITE(*,*) ' % sup                            ', zperx
    WRITE(*,*)

    !
    ! 8. Write optimum in a file
    ! --------------------------
    !
    WRITE(cdum,'(a,1h-,i3.3,1hx,i3.3,1h_,i3.3)') TRIM(covdta),iii,iij,iii*iij-iptx
    OPEN (20,file=cdum)
    WRITE(20,'(a,i5)')'#',iii*iij -iptx
    DO jjc=1,iii*iij-iptx
       WRITE(20,'(a,i5)')'#',jjc
       WRITE(20,'(2i5)')nldi(jjc),nldj(jjc)
       WRITE(20,'(2i5)')nlei(jjc),nldj(jjc)
       WRITE(20,'(2i5)')nlei(jjc),nlej(jjc)
       WRITE(20,'(2i5)')nldi(jjc),nlej(jjc)
       WRITE(20,'(2i5)')nldi(jjc),nldj(jjc)
       WRITE(20,'(2i5)') 9999,9999
       ! Write warnings on standard entry if very few water points (arbitrary <20)
       ! We could test the ratio instead
       IF (icount(jjc).LT.20) THEN
          WRITE(iumout,*)' proc ji=',jjc,' water points:', icount(jjc)
          WRITE(iumout,*) ' ji from ',nldi(jjc), ' to :',nlei(jjc)
          WRITE(iumout,*) ' jj /  mask value for all ji'
          DO jj=nldj(jjc),nlej(jjc)
             WRITE(iumout,900) jj,(INT(zmask(ji,jj)),ji=nldi(jjc),nlei(jjc))
          ENDDO
  900           FORMAT(1x,i4,1x,9(10i1,1x))
       ENDIF
    END DO
    WRITE(20,'(a,i5)')'# vides:',iptx
    DO jjc=1,iptx
       WRITE(20,'(a,i5)')'# vide ',jjc
       WRITE(20,'(2i5)')nldiv(jjc),nldjv(jjc)
       WRITE(20,'(2i5)')nleiv(jjc),nldjv(jjc)
       WRITE(20,'(2i5)')nleiv(jjc),nlejv(jjc)
       WRITE(20,'(2i5)')nldiv(jjc),nlejv(jjc)
       WRITE(20,'(2i5)')nldiv(jjc),nldjv(jjc)
       WRITE(20,'(2i5)')nleiv(jjc),nlejv(jjc)
       WRITE(20,'(2i5)')nldiv(jjc),nlejv(jjc)
       WRITE(20,'(2i5)')nleiv(jjc),nldjv(jjc)
       WRITE(20,'(2i5)') 9999,9999
    END DO
    CLOSE(20)
    IF(iumout .NE. 6) CLOSE(iumout)
    !
    STOP
    !
END PROGRAM mpp_optimiz_nc
