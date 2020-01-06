! Save soil temperatures
       IF (callsoil) THEN
           wtsoil(:,:)=tsoil(:,:)
       ENDIF
! Save fluxrad
       wfluxrad(:)=fluxrad(:)
! Save co2ice if changed
       wco2ice(:)=co2ice(:)
! Save q2
       wq2(:,:)=q2(:,:)
! Save qsurf
       wqsurf(:,:)=qsurf(:,:)
! Save wstar
       wwstar(:)=wstar(:)
