!
!  Apply tracers point sources to the horizontal advection terms,
!  if any, including ICEPLUME type tracers.
!
          IF (LuvSrc(ng)) THEN
            DO is=1,Nsrc(ng)
              IF (INT(PLUME(ng)%dir(is)).ne.0) THEN
                Isrc=SOURCES(ng)%Isrc(is)
                Jsrc=SOURCES(ng)%Jsrc(is)
                IF (((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.          &
     &              ((Jstr.le.Jsrc).and.(Jsrc.le.Jend+1))) THEN
                  IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
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
                    ELSE
                      FX(Isrc,Jsrc)=0.0_r8
                    END IF
                  ELSE
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
                    ELSE
                       FE(Isrc,Jsrc)=0.0_r8
                    END IF
                  END IF
                END IF
              END IF
            END DO
          END IF
