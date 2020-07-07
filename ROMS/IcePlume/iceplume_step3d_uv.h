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
                DO k=1,N(ng)
                  cff1=1.0_r8/(on_u(i,j)*                               &
     &                         0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+     &
     &                                 z_w(i  ,j,k)-z_w(i  ,j,k-1)))
                  u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
                  u(i,j,k,nnew)=u(i,j,k,nnew)+PLUME(ng)%dir(is)*        &
     &              (PLUME(ng)%det(is,k)+PLUME(ng)%ent(is,k)+           &
     &               PLUME(ng)%mB(is,k))*cff1
                END DO
              ELSE
                DO k=1,N(ng)
                  cff1=1.0_r8/(om_v(i,j)*                               &
     &                         0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+     &
     &                                 z_w(i,j  ,k)-z_w(i,j  ,k-1)))
                  v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
                  v(i,j,k,nnew)=v(i,j,k,nnew)+PLUME(ng)%dir(is)*        &
     &              (PLUME(ng)%det(is,k)+PLUME(ng)%ent(is,k)+           &
     &               PLUME(ng)%mB(is,k))*cff1
                END DO
              END IF
            END IF
          END IF
        END DO
      END IF
