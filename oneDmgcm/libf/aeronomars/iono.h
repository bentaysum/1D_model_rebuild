      common/eqF/ o1d_eq,ho2_eq,oh_eq,h_eq,n2d_eq,no2_eq,
     $ o3_eq,no_eq,
     $ cplus_eq,coplus_eq,oplus_eq,n2plus_eq,hplus_eq,co2plus_eq,
     $ o2plus_eq,noplus_eq,nplus_eq,hco2plus_eq

	character*1 o1d_eq(nlayermx),ho2_eq(nlayermx),oh_eq(nlayermx)
      character*1 h_eq(nlayermx),n2d_eq(nlayermx),no2_eq(nlayermx)
      character*1 o3_eq(nlayermx),no_eq(nlayermx),cplus_eq(nlayermx)
      character*1 coplus_eq(nlayermx),oplus_eq(nlayermx)
      character*1 n2plus_eq(nlayermx)
      character*1 hplus_eq(nlayermx),co2plus_eq(nlayermx)
      character*1 o2plus_eq(nlayermx)
      character*1 noplus_eq(nlayermx), nplus_eq(nlayermx) 
      character*1 hco2plus_eq(nlayermx)      

      common/taus/ tauco2,tauo2,tauo3p,tauco,tauh,tauoh,tauho2,tauh2,
     $ tauh2o,tauo1d,tauh2o2,tauo3,taun,tauno,taun2,taun2d,tauno2,
     $ tauco2plus,tauoplus,tauo2plus,taucoplus,taucplus,taunplus,
     $ taunoplus,taun2plus,tauhplus,tauhco2plus

      real*8    tauco2(nreact,nlayermx)
      real*8    tauo2(nreact,nlayermx)
      real*8    tauo3p(nreact,nlayermx)
      real*8    tauco(nreact,nlayermx)
      real*8    tauh(nreact,nlayermx)
      real*8    tauoh(nreact,nlayermx)
      real*8    tauho2(nreact,nlayermx)
      real*8    tauh2(nreact,nlayermx)
      real*8    tauh2o(nreact,nlayermx)
      real*8    tauo1d(nreact,nlayermx)
      real*8    tauh2o2(nreact,nlayermx)
      real*8    tauo3(nreact,nlayermx)
      real*8    taun(nreact,nlayermx)
      real*8    tauno(nreact,nlayermx)
      real*8    taun2(nreact,nlayermx)
      real*8    taun2d(nreact,nlayermx)
      real*8    tauno2(nreact,nlayermx)
      real*8    tauco2plus(nreact,nlayermx)
      real*8    tauoplus(nreact,nlayermx)
      real*8    tauo2plus(nreact,nlayermx)
      real*8    taucoplus(nreact,nlayermx)
      real*8    taucplus(nreact,nlayermx)
      real*8    taunplus(nreact,nlayermx) 
      real*8    taunoplus(nreact,nlayermx)
      real*8    taun2plus(nreact,nlayermx)
      real*8    tauhplus(nreact,nlayermx)
      real*8    tauhco2plus(nreact,nlayermx)
