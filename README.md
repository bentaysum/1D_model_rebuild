# Cleaned 1-D Photochemistry Submodule 

This git repository is initialised with the original model code as recieved from the 
OU, and is altered to provide tracer support for organic chemistry routines involving
CH4 and C2H6. 

## Complaints

Git is arguably the worst software I've ever encountered. This is not the README that
all logical thought processes dictate should be here, but git doesn't work logically,
it operates in some indecipherable way that's either mentally 1000's of years more
advanced than current human thought, or 10000's of years backwards in time in the days
of cavemen. This might be the 100th attempt I've made at creating a flowing TLM branch.
This software is garbage to work with. It's added skip-fulls of stress onto my PhD, 
and the people who make it should truly consider making it at least semi-operatable for
people who haven't got Masters degrees in computing science. Here, I will attempt to
record the building of the TLM. It won't work; it will become messy, convoluted, 
non-navigatable, and at this point all I can see is: I don't frankly care anymore. Thank
you github, for showing me how inept I truly am at the should-be-easy process of version
control. I'm earnestly trying my best here.


## 03/04/2020

Attempting to make the TLM operate over the original number of discrete time-steps (48 physical)
to enable the photochemistry routine to function over 3 chemistry time-steps, and also the 5 sub-
timesteps within each chemistry timestep. 

## 08/04/2020

Created a quick routine that computes the perturbation vector within the 1-D Model

	- Lines 2358 in physiq.F
	- pertvector in TLMvars and TLM_initialise
	- Line 1344 in physiq.F

Bug stomped. Routine operation. 

# VALIDATION PENDING

## 07/04/2020

dccn_dpq and dcc0_dpq are no longer global save variables.
dHOX_dpq and dHOX0_dpq are no longer global saves.
dOX_dpq and dOX0_dpq are no longer global saves.

Linearised variables initialised outside of chemistry loop.

Input into photochemistry matrix is outside of the chemistry loop.

Operational (no errors) with nsteps = 48 per sol!

## 31/03/2020

I have inserted the files:

#### libf/aeronomars/tlm_ox.F90

Linearised photochemistry calculations concerning the advancement of odd-oxygen (O,O3,O1D)
compounds during daylight conditions. Calculates the advancement of CC[OX].

#### libf/aeronomars/tlm_hox.F90

Linearised photochemistry calculations concerning the advancement of odd-hydrogen (H,OH,HO2)
compounds at all model iterations. Advancement of CC[HOx] is *not* handled here.

#### libf/aeronomars/tlm_sibem.F90

Covers the advancement of all other chemical species. Advancement of CC[HOX] *is* handled
here.

#### libf/phymars/TLMvars.F90

A module handling the global variables:

- t_X [nqmx*nlayermx] : ordering of the chemical species within TLM vector space
- tlm [nqmx*nlayermx,nqmx*nlayermx] : the Tangent Linear Model matrix
- tlm_photo [nqmx*nlatermx,nqmx*nlayermx]: photochemical component of the tangent linear matrix
- dhox_dpq [nlayermx,nqmx*nlayermx] : linearised advancement of HOx number denisty
- dox_dpq [nlayermx,nqmx*nlayermx] : linearised advancement of Ox number density
- dccn_dpq [nqmx*nlayermx,nqmx*nlayermx] : linearised advancement of the CC (number density) vector
- dcc0_dpq [nqmx*nlayermx,nqmx*nlayermx] : linearised initial number density entering routine
- dhox0_dpq [nqmx*nlayermx,nqmx*nlayermx] : " " HOx number density " "
- dox0_dpq [nqmx*nlayermx,nqmx*nlayermx] : " " OX number density " "
- TLM_ident [nqmx*nlayermx,nqmx*nlayermx] : Identity matrix of the TLM
- Avmr [nlayermx,nqmx] : mass to volume mixing ratio conversion factors

#### libf/phymars/TLM_vdif.F90

Linearised vertical transport routine.

#### libf/phymars/TLM_initialise.F

Initialises the tangent linear model variables via allocation statements.

#### libf/phymars/TLM_ADJ_out.F90

Stores each model iterations TLM matrix into the tangent_matrix [nqmx*nlayermx,nqmx*nlayermx,ndt]
(for ndt model time-steps) to be saved as a binary output file.

## 01/04/2020 

Routine now operates with a bash script that cycles through a user defined number of scalar 
perturbations and stores them in the university datastore. These 1-D model outputs can be used to 
assess the validity of the Tangent Linear Model code. The cycles perturb the species at the same
time index that the TLM is supposed to be activated on.

To access TLM validity, set up the bash script to specify the:
	- TRACER [tracer to pertrub]
	- LOW [smallest magnitude of perturbation]
	- HIGH [largest magnitude of perturbation]
	- L_1 [ lowest layer in model for perturbation]
	- L_2 [ upper layer in model for perturbation]

and also the day0 value. Then, activate the call_tlm callkey in callphys.def, DEACTIVATE THE 
PERTLOOP KEY, and you have yourself:
	- A tangent linear matrix (TLM.bin)
	- Integrated model perturbation files
	- A control run 
