      INTEGER wday_ini
      REAL wtsurf(ngridmx)  ! input only ay firstcall - output
      REAL wtsoil(ngridmx,nsoilmx)
      REAL wisoil(ngridmx,nsoilmx)  !! new soil scheme
      REAL wdsoil(ngridmx,nsoilmx)   !! new soil scheme
      REAL wco2ice(ngridmx)
      REAL wemis(ngridmx)
      REAL wqsurf(ngridmx,nqmx)
      REAL wq2(ngridmx,nlayermx+1)
      REAL wwstar(ngridmx)
      REAL wfluxrad(ngridmx)
      REAL output_tab2d(ngridmx,n2d)
      REAL output_tab3d(ngridmx,nlayer,n3d)
      REAL sensheat(ngridmx)    !! pour LES avec isfflx!=0
      !REAL ustar(ngridmx)    !! pour LES avec isfflx!=0 !! already defined
      LOGICAL flag_LES     !! pour LES avec isfflx!=0
      REAL qsurfice(ngridmx) !! pour diagnostics
      real alpha,lay1 ! coefficients for building layers
      integer iloop
      INTEGER tracerset    !!! this corresponds to config%mars
      CHARACTER (len=20) :: wtnom(nqmx) ! tracer name

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JF 
      REAL qsurfdust(ngridmx)         ! useful for dust diagnostics
      REAL TAU_lay(ngridmx)           ! true opacity (it's not e reference life tauref)
      REAL dsodust(ngridmx,nlayermx)
      REAL zdqsdif_diag(ngridmx)      ! useful for lifting diagnostics                               
      REAL zdqssed_diag(ngridmx)      ! useful for sedimentation diagnostics                         
      REAL pdq_diag(ngridmx)          ! useful for dust perturbation diagnostics   
      REAL dustot(ngridmx)            ! Total mass of dust integrated along vertical axe (kg/m2)
      REAL zdqnorm(ngridmx,nlayermx,2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! JF
