# Organic 1-D Photochemistry Submodule 

This git repository is initialised with the original model code as recieved from the 
OU, and is altered to provide tracer support for organic chemistry routines involving
CH4 and C2H6.

Model version used in production of data used in the AGU submitted publication:
"Photochemistry of Methane and Ethane in the Martian Atmosphere".

## Compilation and successful runs
Model has been tailored to work with the gfortran compiler via the command inside of the oneDmgcm directory:

`./makegcm_gfortran -d 1x1x25 -p mars testphys1d`

To enable successful compilation:

### oneDmgcm/makegcm_gfortran

Specify locations of netCDF libraries through the NCDFLIB and NCDFINC variables.

Number of tracers can be selected via the -t XX flag at compile statement where XX is an integer (27 for only methane oxidation, 57 for methane and ethane oxidation).

### oneDmgcm/traceur.def

Definition of tracer names. Ensure the integer in the first row == XX in the compile input. Do not alter tracer names.

### oneDmgcm/callphys.def

datadir : define location of the mgcm-datafile directory.

### oneDmgcm/libf/mcd_trac.F

edit FILE_NAME variable to specify location of the mcd-datafile directory.

### oneDmgcm/libf/mcd_phys.F

edit FILE_NAME variable to specify location of the mcd-datafile directory. 

## Operating the 1-D Organic model

Model operation is chiefly performed through the run.def and callphys.def input files located within the oneDmgcm directory.

### oneDmgcm/callphys.def

#### General 
- diurnal : controls solar zenith calculations; .TRUE. for all experiments.

- seasonal : controls calculations for solar longitude calculation; .TRUE. for annual runs, .FALSE. for shorter temporal studies.

- save_rates : controls output - .TRUE. for rate coefficients, .FALSE. for mass mixing ratios (mmr)

- mcd_profiles : controls tracer profile interpolation routines from MCDv5.3 dataset; .TRUE. for all experiments.

- long_mean : interpolates from a longitudinal mean MCDv5.3 grid; .TRUE. for all experiments (issues with time-keeping+solar zenith angle calculations when .FALSE.)

- mcd_dayfin : number of sols for the model to run for prior to turning of the tracer PROFILE interpolation schemes; > ndt (number of sols whole model) for all experiments.

- mcd_lt : local time to shut profile interpolation down at (arbitrary in this study, must be 0.0 <= mcd_lt < 24.0 however)

#### Organic releases
- plume : instantaneously deposit an organic species; .False. for annual runs, .True. for shorter temporal studies.
- plume_day : number of days for model to run for prior to plume release, also shutting down MCD tracer interpolations; 10.
- plume_lt : local time to release organic tracers at; 0.5.
- p_low : lowest level of insertion; 1
- p_high : highest layer of insertion; 25

### oneDmgcm/run.def
- day0 : a table mapping integer values to solar longitude values is found in the user_manual.pdf, -10 for an annual run.
- time : 0 for all studies
- ndt : number of sols to run for; 678 for an annual run, 20 for a shorter temporal study.
- latitude : model latitude; 30., 2.5, and -30. in our works.




