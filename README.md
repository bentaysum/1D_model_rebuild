# Cleaned 1-D Photochemistry Submodule 

This git repository is initialised with the original model code as recieved from the 
OU, and is altered to provide tracer support for organic chemistry routines involving
CH4 and C2H6. 

## 07/01/2020

Module now interpolates from the MCD v5.3 offline look-up tables continuously for
profiles of CO2, CO, O2, H2 and H2O vapour as mixing ratios, and profiles of the
physical variables Temperature, Q2, u, v and surface pressure. 

Module now equipped with the callkeys mcd_dayfin and mcd_ltfin to indicate the 
number of days the tracer interpolation routine runs for, and what local time (lt)
we want the module to stop interpolating these tracer profiles at.

tracer_release.F introduced to read concentration and tracer parameters from the plume.def
file that will tell the model what gas will be inserted across layers p_low to p_high [defined
in callphys.def] and at what VMR. 

## ONGOING 
finish the tracer_release.F file. Make it work for PROFILE concentrations, and also
COLUMN VMR abundances concentrated in between layers p_low and p_high
