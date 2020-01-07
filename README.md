# Cleaned 1-D Photochemistry Submodule 

This git repository is initialised with the original model code as recieved from the 
OU, and is altered to provide tracer support for organic chemistry routines involving
CH4 and C2H6. 

## 07/01/2020

Module now interpolates from the MCD v5.3 offline look-up tables continuously for
profiles of CO2, CO, O2, H2 and H2O vapour as mixing ratios, and profiles of the
physical variables Temperature, Q2, u, v and surface pressure. 

Next step: establish a cut-off routine for tracer profile interpolation to enable
the study of the natural evolution of the tracers. 
