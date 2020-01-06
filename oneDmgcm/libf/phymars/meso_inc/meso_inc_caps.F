      !!! This has to go after initracer which change dryness and watercaptag

      !!!!!! MARS MESOSCALE MODELING
      !!!!!! TEST TEST TEST TEST  AS+JBM 28/02/11
      !!!!!! TEST TEST TEST TEST  AS+JBM 28/02/11
      !!!!!! TEST TEST TEST TEST  AS+JBM 28/02/11
      !!!! 
      !!!! principe: une option 'caps=T' specifique au mesoscale
      !!!! ... en vue d'un meso_initracer ????
      !!!! 
      !!!! depots permanents => albedo TES du PDS
      !!!! depots saisonniers => alb_surfice (~0.4, cf plus bas)
      !!!!     [!!!! y compris pour les depots saisonniers sur les depots permanents]
      !!!!
      !!!! --> todo: il faut garder les depots saisonniers qui viennent
      !!!!           du GCM lorsqu'ils sont consequents
      !!!! 
      !IF ( caps .and. (igcm_h2o_ice .ne. 0) ) THEN

      IF ( caps .and. water ) THEN
          PRINT *, 'OVERWRITING watercaptag DEFINITION in INITRACER'
          PRINT *, 'lat > lat_lim et alb > alb_lim => watercaptag=T'
          PRINT *, 'ind for water ice: ', igcm_h2o_ice 
          !! Perennial H20 north cap defined by watercaptag=true (allows surface to be
          !! hollowed by sublimation in vdifc).
          do ig=1,ngridmx
            qsurf(ig,igcm_h2o_ice)=0.  !! on jette les inputs GCM
            if ( ( lati(ig)*180./pi .gt. 70. ) .and.
     .           ( albedodat(ig) .ge. 0.26   ) )  then
                    watercaptag(ig)  = .true.
                    dryness(ig)      = 1.
            else
                    watercaptag(ig)  = .false.
                    dryness(ig)      = 1.
            endif  ! (lati, albedodat)
          end do ! (ngridmx)
      ELSE  ! (caps)
          print *,'Blork !!!'
          print *,'caps=T avec water=F ????'
      ENDIF ! (caps)
