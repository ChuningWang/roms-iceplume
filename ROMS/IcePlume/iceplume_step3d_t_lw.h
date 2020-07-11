!
!  Apply mass point sources (volume vertical influx), if any, including
!  iceplume type tracers.
!
          IF (LwSrc(ng)) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng)) THEN
                IF ((Vadvection(itrc,ng)%MPDATA).or.                    &
     &              (Vadvection(itrc,ng)%HSIMT)) THEN
                  LapplySrc=(IstrUm2.le.Isrc).and.                      &
     &                      (Isrc.le.Iendp2i).and.(j.eq.Jsrc)
                ELSE
                  LapplySrc=(IstrR.le.Isrc).and.                        &
     &                      (Isrc.le.IendR).and.(j.eq.Jsrc)
                END IF
                IF (LapplySrc) THEN
                  cff=0.0_r8
                  DO k=1,N(ng)
                    cff=cff+SOURCES(ng)%Qsrc(is,k)*                     &
     &                      SOURCES(ng)%Tsrc(is,k,itrc)
                    cff=cff+PLUME(ng)%det   (is,k     )*                &
# ifdef ICEPLUME_DET_NEUTRAL
     &                      PLUME(ng)%detTrc(is,k,itrc)+                &
# else
     &                      PLUME(ng)%trc   (is,  itrc)+                &
# endif
     &                      PLUME(ng)%ent   (is,k     )*                &
     &                      PLUME(ng)%trcAm (is,k,itrc)+                &
     &                      PLUME(ng)%mB    (is,k     )*                &
     &                      PLUME(ng)%trcB  (is,  itrc)
# ifdef SPLINES_VDIFF
                    IF (.not.((Hadvection(itrc,ng)%MPDATA).and.         &
     &                 (Vadvection(itrc,ng)%MPDATA))) THEN
                      cff=cff*oHz(Isrc,Jsrc,k)
                    END IF
# endif
                    FC(Isrc,k)=FC(Isrc,k)-cff
                  END DO
                END IF
              END IF
            END DO
          END IF
