      IF (flag_LES) THEN
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdif(ig,l,iq)
              ENDDO
            ENDDO
           ENDDO
           DO iq=1, nq
              DO ig=1,ngrid
                 dqsurf(ig,iq)=dqsurf(ig,iq) + zdqsdif(ig,iq)
              ENDDO
           ENDDO
      ENDIF
