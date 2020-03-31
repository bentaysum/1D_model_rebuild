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

## 08/01/202

tracer_release.F finished for profile perturbations. Reads plume.def and deposits the specified
abundance of the listed tracers instantaneously and homogeneously through layer p_low to p_high,
provided in the callphys.def file. Once routine is called, the module immediately stops interpolating
tracer profiles from the MCD v5.3 offline look-up table.
