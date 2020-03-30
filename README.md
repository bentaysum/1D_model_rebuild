# Cleaned 1-D Photochemistry Submodule 

This git repository is initialised with the original model code as recieved from the 
OU, and is altered to provide tracer support for organic chemistry routines involving
CH4 and C2H6. 

## Version for AGU Submission (Photochemistry of Methane and Ethane)

We compile our 1-D model code using the gfortran compiler, through the command:

./makegcm_gfortran -d 1x1x25 -p mars testphys1d

inside of the /oneDmgcm directory to produce the executable testphys1d.e .

## netCDF

Insert location of netcdf libraries and include files in the makegcm_gfortran script:

setenv NCDFLIB /usr/lib64/gfortran/modules
setenv NCDFINC /usr/include 

alter these lines as required.
 
## Atmospheric parameters

1-D Model uses /oneDmgcm/libf/phymars/mcd_phys.F to linearly interpolate atmospheric
parameters with local time, solar longitude, altitude, and latitude.

The variables FILE_NAME must be changed to the relevant name pointing to the mcd_physicals.nc
dataset stored on your machine. The iaervar == 25 conditional is a build in progress and 
ignorable for this work.

## Tracer volume mixing ratios

1-D model uses /oneDmgcm/libf/phymars/mcd_trac.F to similarly interpolate tracer vmrs of 
CO2, CO, O2, H2O vapour and H2 from the Mars Climate Database using a mean molecular mass
of air of 43.34.

Similarly, specify the FILE_NAME variable to point to the location of the mcd_tracervmrs.nc
dataset on your machine.

## callphys.def 

We specify relevant callkeys in callphys.def that are changed to produce our data/get the model
working on different machines. All keys not mentioned should be unchanged.

The 1-D model makes use of datafiles from the MGCM, located in the mgcm-datafile directory,
specified in callphys.def. Change the datadir variable in here accordingly.

season : set to true for annual runs, off for single/double digit sol runs

chemfluxes : false for mass mixing ratio outputs, true for the reaction rate coefficients.

plume : true deposits a user specified abundance of gas described in plume.def in the model, 
       which shuts down the tracer interpolation routines at the point of activation. p_low, 
       p_high control the altitude range of the release. These don't vary in the paper. plume_lt
       specifies the local time of deposit, and plume_day specifies how many sols must elapse 
       prior to activation. Again, invariant in paper.

## run.def

day0 : initial date of model (perihelion and aphelion locations specified), set to -10 for annual
       run
ndt : 20 for short runs, 679 for annual runs.

latitude : specifies model latitude.

## plume.def
specifies tracer to release, and at what concentration (when plume == .True.)
