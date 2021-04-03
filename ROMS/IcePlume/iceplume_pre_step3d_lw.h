!
!  Apply mass point sources (volume vertical influx), if any, including
!  iceplume type tracers.
!  (Notice the j-loop is pipelined above).
!
          IF (LwSrc(ng)) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
     &            ((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.            &
     &            (j.eq.Jsrc)) THEN
                cff=0.0_r8
                DO k=1,N(ng)
                  cff=cff+SOURCES(ng)%Qsrc(is,k)*                       &
     &                    SOURCES(ng)%Tsrc(is,k,itrc)
                  cff=cff+PLUME(ng)%det   (is,k     )*                  &
#  ifdef ICEPLUME_DET_NEUTRAL
     &                    PLUME(ng)%detTrc(is,k,itrc)+                  &
#  else
     &                    PLUME(ng)%trc   (is,  itrc)+                  &
#  endif
     &                    PLUME(ng)%ent   (is,k     )*                  &
     &                    PLUME(ng)%trcAm (is,k,itrc)+                  &
     &                    PLUME(ng)%mB    (is,k     )*                  &
     &                    PLUME(ng)%trcB  (is,k,itrc)
                  FC(Isrc,k)=FC(Isrc,k)-cff
                END DO
              END IF
            END DO
          END IF
