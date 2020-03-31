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

