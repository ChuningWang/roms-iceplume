!
!  Apply mass point sources (volume vertical influx), if any, including
!  ICEPLUME type tracers.
!
        IF (LwSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            ii=SOURCES(ng)%Isrc(is)
            jj=SOURCES(ng)%Jsrc(is)
            IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                  &
     &          ((JstrR.le.jj).and.(jj.le.JendR)).and.                  &
     &          (j.eq.jj)) THEN
              DO k=1,N(ng)
                W(ii,jj,k)=W(ii,jj,k-1)-                                &
     &               (Huon(ii+1,jj,k)-Huon(ii,jj,k)+                    &
     &                Hvom(ii,jj+1,k)-Hvom(ii,jj,k))+                   &
     &                SOURCES(ng)%Qsrc(is,k)
                W(ii,jj,k)=W(ii,jj,k)+                                  &
     &                     (PLUME(ng)%det(is,k)+PLUME(ng)%ent(is,k)+    &
     &                      PLUME(ng)%mB(is,k))
              END DO
            END IF
          END DO
        END IF
