!
!  Apply tracers point sources to the horizontal advection terms,
!  if any, including iceplume type tracers.
!
          IF (LuvSrc(ng)) THEN
            DO is=1,Nsrc(ng)
              IF (INT(PLUME(ng)%dir(is)).ne.0) THEN
                Isrc=SOURCES(ng)%Isrc(is)
                Jsrc=SOURCES(ng)%Jsrc(is)
                IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                  IF ((Hadvection(itrc,ng)%MPDATA).or.                  &
     &                (Hadvection(itrc,ng)%HSIMT)) THEN
                    LapplySrc=(IstrUm2.le.Isrc).and.                    &
     &                        (Isrc.le.Iendp3).and.                     &
     &                        (JstrVm2.le.Jsrc).and.                    &
     &                        (Jsrc.le.Jendp2i)
                  ELSE
                    LapplySrc=(Istr.le.Isrc).and.                       &
     &                        (Isrc.le.Iend+1).and.                     &
     &                        (Jstr.le.Jsrc).and.                       &
     &                        (Jsrc.le.Jend)
                  END IF
                  IF (LapplySrc) THEN
                    IF (LtracerSrc(itrc,ng)) THEN
                      FX(Isrc,Jsrc)=SOURCES(ng)%Qsrc(is,k     )*        &
     &                              SOURCES(ng)%Tsrc(is,k,itrc)
                      FX(Isrc,Jsrc)=FX(Isrc,Jsrc)+PLUME(ng)%dir(is)*(   &
# ifdef ICEPLUME_DET_NEUTRAL
     &                 PLUME(ng)%det(is,k)*PLUME(ng)%detTrc(is,k,itrc)+ &
# else
     &                 PLUME(ng)%det(is,k)*PLUME(ng)%trc   (is,  itrc)+ &
# endif
     &                 PLUME(ng)%ent(is,k)*PLUME(ng)%trcAm (is,k,itrc)+ &
     &                 PLUME(ng)%mB (is,k)*PLUME(ng)%trcB  (is,  itrc))
# ifdef MASKING
                    ELSE
                      IF ((rmask(Isrc  ,Jsrc).eq.0.0_r8).and.           &
     &                    (rmask(Isrc-1,Jsrc).eq.1.0_r8)) THEN
                        FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                &
     &                                t(Isrc-1,Jsrc,k,3,itrc)
                      ELSE IF ((rmask(Isrc  ,Jsrc).eq.1.0_r8).and.      &
     &                         (rmask(Isrc-1,Jsrc).eq.0.0_r8)) THEN
                        FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                &
     &                                t(Isrc  ,Jsrc,k,3,itrc)
                      END IF
# endif
                    END IF
                  END IF
                ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
                  IF ((Hadvection(itrc,ng)%MPDATA).or.                  &
     &                (Hadvection(itrc,ng)%HSIMT)) THEN
                    LapplySrc=(IstrUm2.le.Isrc).and.                    &
     &                        (Isrc.le.Iendp2i).and.                    &
     &                        (JstrVm2.le.Jsrc).and.                    &
     &                        (Jsrc.le.Jendp3)
                  ELSE
                    LapplySrc=(Istr.le.Isrc).and.                       &
     &                        (Isrc.le.Iend).and.                       &
     &                        (Jstr.le.Jsrc).and.                       &
     &                        (Jsrc.le.Jend+1)
                  END IF
                  IF (LapplySrc) THEN
                    IF (LtracerSrc(itrc,ng)) THEN
                      FE(Isrc,Jsrc)=SOURCES(ng)%Qsrc(is,k     )*        &
     &                              SOURCES(ng)%Tsrc(is,k,itrc)
                      FE(Isrc,Jsrc)=FE(Isrc,Jsrc)+PLUME(ng)%dir(is)*(   &
# ifdef ICEPLUME_DET_NEUTRAL
     &                 PLUME(ng)%det(is,k)*PLUME(ng)%detTrc(is,k,itrc)+ &
# else
     &                 PLUME(ng)%det(is,k)*PLUME(ng)%trc   (is,  itrc)+ &
# endif
     &                 PLUME(ng)%ent(is,k)*PLUME(ng)%trcAm (is,k,itrc)+ &
     &                 PLUME(ng)%mB (is,k)*PLUME(ng)%trcB  (is,  itrc))
# ifdef MASKING
                    ELSE
                      IF ((rmask(Isrc,Jsrc  ).eq.0.0_r8).and.           &
     &                    (rmask(Isrc,Jsrc-1).eq.1.0_r8)) THEN
                        FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                &
     &                                t(Isrc,Jsrc-1,k,3,itrc)
                      ELSE IF ((rmask(Isrc,Jsrc  ).eq.1.0_r8).and.      &
     &                         (rmask(Isrc,Jsrc-1).eq.0.0_r8)) THEN
                        FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                &
     &                                t(Isrc,Jsrc  ,k,3,itrc)
                      END IF
# endif
                    END IF
                  END IF
                END IF
              END IF
            END DO
          END IF
