      SUBROUTINE tracer_release(dq)

      IMPLICIT NONE

#include "dimensions.h"
#include "dimphys.h"
#include "comgeomfi.h"
#include "surfdat.h"
#include "comsoil.h"
#include "comdiurn.h"
#include "callkeys.h"
#include "comcstfi.h"
#include "planete.h"
#include "comsaison.h"
#include "control.h"
#include "dimradmars.h"
#include "comg1d.h"
#include "tracer.h"
#include "nlteparams.h"

#include "chimiedata.h"
#include "param.h"
#include "param_v4.h"
#include "conc.h"

!     INPUT
      real dq(nlayermx,nqmx) ! Tracer tendency from 1-D model

!     LOCAL       
      real dq_plume(nlayermx,nqmx) ! Tracer tendencies to be added due to plume 
      
      character(len=*), parameter :: filename = "plume.def" 
      
      
      OPEN(unit=100,file = filename)
      WRITE(*,*) "OH"
      STOP
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      END 