      COMMON/surfdat/albedodat(ngridmx),                                &
     &   phisfi(ngridmx),albedice(2),emisice(2),emissiv,                &
     &   TESice_Ncoef,TESice_Scoef,                                     &
     &   iceradius(2) , dtemisice(2),                                   &
     &   zmea(ngridmx),zstd(ngridmx),                                   &
     &   zsig(ngridmx),zgam(ngridmx),zthe(ngridmx),                     &
     &   z0(ngridmx),z0_default, albedo_h2o_ice, inert_h2o_ice,         &
     &   frost_albedo_threshold

      COMMON/surfdatl/TESicealbedo,watercaptag,temptag


      real albedodat ! albedo of bare ground
      real phisfi ! geopotential at ground level
      real albedice ! default albedo for ice (1: North H. 2: South H.)
      real emisice ! ice emissivity; 1:Northern hemisphere 2:Southern hemisphere
      real emissiv ! emissivity of bare ground
      logical TESicealbedo ! use TES ice cap albedoes (if set to .true.)
      logical watercaptag(ngridmx) ! flag for water ice surface
      
      logical temptag !temp tag for water caps
      
      real albedo_h2o_ice ! water ice albedo
      real inert_h2o_ice ! water ice thermal inertia
      real frost_albedo_threshold ! water frost thickness on the ground (kg.m^-2, ie mm)
      real TESice_Ncoef ! coefficient for TES ice albedo in Northern hemisphere
      real TESice_Scoef ! coefficient for TES ice albedo in Southern hemisphere
      real iceradius , dtemisice
      real zmea,zstd,zsig,zgam,zthe
      real z0 ! surface roughness length (m)
      real z0_default ! default (constant over planet) surface roughness (m)
