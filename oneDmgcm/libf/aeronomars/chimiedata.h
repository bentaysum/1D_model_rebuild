! -----------------------------------------
! ---  Data for photolysis lookup table ---
! -----------------------------------------

      integer nd, nozo, nr, nsza, ntemp, ntau

      parameter (nd    = 29)
      parameter (nozo  = 7)
      parameter (nr    = nd + 28)
      parameter (nsza  = 27)
      parameter (ntemp = 4)
      parameter (ntau  = 8)

      real kb
      parameter (kb = 1.3806e-23)

      common/chimiedata/jphot,colairtab,table_ozo

      real jphot(ntemp,nsza,0:200,nozo,ntau,nd)
      real colairtab(0:200)
      real szatab(nsza)
      real table_ozo(nozo)
      real tautab(ntau)

      data szatab/0.,  5., 10., 15., 20., 25.,                          &
     &            30., 35., 40., 45., 50., 55.,                         &
     &            60., 65., 70., 75., 80., 82.,                         &
     &            84., 86., 88., 90., 91., 92.,                         &
     &            93., 94., 95./

      data tautab/0., 0.2, 0.4, 0.6, 0.8, 1., 2., 4./
