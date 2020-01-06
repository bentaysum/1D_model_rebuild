c-----------------------------------------------------------------------
c INCLUDE 'sponge.h'

      COMMON/com_sponge_l/callsponge
      common/com_sponge_i/mode_sponge,nsponge
      common/com_sponge_r/tetasponge

      LOGICAL   callsponge  ! do we use a sponge on upper layers
      INTEGER   mode_sponge ! sponge mode
      INTEGER nsponge ! number of sponge layers 
      REAL  tetasponge  ! sponge time scale (s) at topmost layer
c-----------------------------------------------------------------------
