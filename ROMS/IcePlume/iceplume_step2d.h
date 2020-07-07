!
!-----------------------------------------------------------------------
!  Apply momentum transport point sources (like river runoff), if any,
!  including ICEPLUME type tracers.
!-----------------------------------------------------------------------
!
      IF (LuvSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          IF (INT(PLUME(ng)%dir(is)).ne.0) THEN
            i=SOURCES(ng)%Isrc(is)
            j=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.i).and.(i.le.IendR)).and.                    &
     &          ((JstrR.le.j).and.(j.le.JendR))) THEN
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
                cff=1.0_r8/(on_u(i,j)*                                  &
     &                      0.5_r8*(zeta(i-1,j,knew)+h(i-1,j)+          &
     &                              zeta(i  ,j,knew)+h(i  ,j)))
                ubar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
                ubar(i,j,knew)=ubar(i,j,knew)+                          &
     &            PLUME(ng)%trs(is)*PLUME(ng)%dir(is)*cff
              ELSE
                cff=1.0_r8/(om_v(i,j)*                                  &
     &                      0.5_r8*(zeta(i,j-1,knew)+h(i,j-1)+          &
     &                              zeta(i,j  ,knew)+h(i,j  )))
                vbar(i,j,knew)=SOURCES(ng)%Qbar(is)*cff
                vbar(i,j,knew)=vbar(i,j,knew)+                          &
     &            PLUME(ng)%trs(is)*PLUME(ng)%dir(is)*cff
              END IF
            END IF
          END IF
        END DO
      END IF
