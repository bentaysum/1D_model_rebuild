#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Definitions de Macros pour Make
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# create_make_gcm arguments: dyn3d grid bibio aeronomars phymars

# Repertoires :

GCM     = /home/s1215319/mgcm/oneDmgcm
MACHINE = sun4
LIBF    = $(GCM)/libf
LIBO    = $(LIBOGCM)/$(MACHINE)
LOCAL_DIR=/home/s1215319/mgcm/oneDmgcm
BIBIO    = $(LIBF)/bibio
AERONOMARS    = $(LIBF)/aeronomars
FILTRE   = filtre
PHYS  = 
DYN  = $(DYNTYPE)
LIBPHY = $(LIBO)/libphy$(PHYS).a
DIRMAIN= $(DYN)$(DIM)d
RM=rm

OPLINK = 

# Les differentes librairies pour l"edition des liens:

dyn3d      = $(LIBO)/libdyn3d.a $(LIBO)/lib$(FILTRE).a
dyn2d      = $(LIBO)/libdyn2d.a
dyn1d      = $(LIBO)/libdyn1d.a
L_DYN      = -l$(DYN)$(DIM)d
L_FILTRE   = -l$(FILTRE)
L_PHY = -lphy$(PHYS) 
L_BIBIO    = -lbibio
L_AERONOMARS    = -laeronomars
L_ADJNT    =

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Option de compilation FORTRAN
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMPILE = $(F77) $(OPTIM) $(INCLUDE) -c
COMPILE90 = $(F90) $(OPTIM90) $(INCLUDE) -c
COMPTRU90 = $(F90) $(OPTIMTRU90) $(INCLUDE) -c
LINK    = 
AR      = 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creation des differents executables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Executables:
# ------------

PROG = code

main : $(DYN) bibio phys aeronomars $(OPTION_DEP) 
	cd $(LIBO) ; $(RANLIB) lib*.a ; cd $(GCM) ;\
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/$(DIRMAIN)/$(SOURCE) -o $(PROG).o ; \
	$(LINK) $(PROG).o -L$(LIBO) $(L_DYN) $(L_ADJNT) $(L_FILTRE) $(L_AERONOMARS) $(L_PHY) $(L_BIBIO) $(L_DYN) $(L_AERONOMARS) $(OPLINK) $(OPTION_LINK) -o $(LOCAL_DIR)/$(PROG).e ; $(RM) $(PROG).o 

dyn : $(LIBO)/libdyn$(DIM)d.a $(FILTRE)$(DIM)d

olddyn : $(LIBO)/libolddyn$(DIM)d.a $(FILTRE)$(DIM)d

phys : $(LIBPHY)

bibio : $(LIBO)/libbibio.a

aeronomars : $(LIBO)/libaeronomars.a

adjnt : $(LIBO)/libadjnt.a

$(FILTRE)3d : $(LIBO)/lib$(FILTRE).a

$(FILTRE)2d :

$(FILTRE)1d :

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Contenu des differentes bibliotheques
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory aeronomars
#=======================================================================

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(concentrations.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(conduction.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(deposition.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(dtridgl.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(hrtherm.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(intrplf.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(inv.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(jthermcalc_e107.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(jthermcalc.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(moldiffcoeff.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(moldiffcoeff_red.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(moldiff.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(molvis.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(paramfoto_compact.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(param_read_e107.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(param_read.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(perosat.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(photochemistry.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(read_phototable.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(surfacearea.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(thermosphere.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(write_ratecoefficient.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(calchim.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(chemthermos.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(euvheat.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(hcl_deposition.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(heterogeneous_chlorine.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(inichim_newstart.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(moldiff_red.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(tgo_profiles.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(TLM_chem.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(tlm_clox.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(tlm_hox.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(tlm_nox.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(tlm_ox.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(tlm_sibem.o)

.PRECIOUS	: $(LIBO)/libaeronomars.a


# Compilation des membres de la bibliotheque libaeronomars.a

$(LIBO)/libaeronomars.a(concentrations.o) : $(LIBF)/aeronomars/concentrations.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a concentrations.o ; \
	$(COMPILE) $(LIBF)/aeronomars/concentrations.F ; \
	$(AR) r $(LIBO)/libaeronomars.a concentrations.o ; $(RM) concentrations.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(conduction.o) : $(LIBF)/aeronomars/conduction.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a conduction.o ; \
	$(COMPILE) $(LIBF)/aeronomars/conduction.F ; \
	$(AR) r $(LIBO)/libaeronomars.a conduction.o ; $(RM) conduction.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(deposition.o) : $(LIBF)/aeronomars/deposition.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a deposition.o ; \
	$(COMPILE) $(LIBF)/aeronomars/deposition.F ; \
	$(AR) r $(LIBO)/libaeronomars.a deposition.o ; $(RM) deposition.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(dtridgl.o) : $(LIBF)/aeronomars/dtridgl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a dtridgl.o ; \
	$(COMPILE) $(LIBF)/aeronomars/dtridgl.F ; \
	$(AR) r $(LIBO)/libaeronomars.a dtridgl.o ; $(RM) dtridgl.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(hrtherm.o) : $(LIBF)/aeronomars/hrtherm.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a hrtherm.o ; \
	$(COMPILE) $(LIBF)/aeronomars/hrtherm.F ; \
	$(AR) r $(LIBO)/libaeronomars.a hrtherm.o ; $(RM) hrtherm.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(intrplf.o) : $(LIBF)/aeronomars/intrplf.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a intrplf.o ; \
	$(COMPILE) $(LIBF)/aeronomars/intrplf.F ; \
	$(AR) r $(LIBO)/libaeronomars.a intrplf.o ; $(RM) intrplf.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(inv.o) : $(LIBF)/aeronomars/inv.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a inv.o ; \
	$(COMPILE) $(LIBF)/aeronomars/inv.F ; \
	$(AR) r $(LIBO)/libaeronomars.a inv.o ; $(RM) inv.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(jthermcalc_e107.o) : $(LIBF)/aeronomars/jthermcalc_e107.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a jthermcalc_e107.o ; \
	$(COMPILE) $(LIBF)/aeronomars/jthermcalc_e107.F ; \
	$(AR) r $(LIBO)/libaeronomars.a jthermcalc_e107.o ; $(RM) jthermcalc_e107.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(jthermcalc.o) : $(LIBF)/aeronomars/jthermcalc.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a jthermcalc.o ; \
	$(COMPILE) $(LIBF)/aeronomars/jthermcalc.F ; \
	$(AR) r $(LIBO)/libaeronomars.a jthermcalc.o ; $(RM) jthermcalc.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(moldiffcoeff.o) : $(LIBF)/aeronomars/moldiffcoeff.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a moldiffcoeff.o ; \
	$(COMPILE) $(LIBF)/aeronomars/moldiffcoeff.F ; \
	$(AR) r $(LIBO)/libaeronomars.a moldiffcoeff.o ; $(RM) moldiffcoeff.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(moldiffcoeff_red.o) : $(LIBF)/aeronomars/moldiffcoeff_red.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h \
$(LIBF)/aeronomars/diffusion.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a moldiffcoeff_red.o ; \
	$(COMPILE) $(LIBF)/aeronomars/moldiffcoeff_red.F ; \
	$(AR) r $(LIBO)/libaeronomars.a moldiffcoeff_red.o ; $(RM) moldiffcoeff_red.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(moldiff.o) : $(LIBF)/aeronomars/moldiff.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a moldiff.o ; \
	$(COMPILE) $(LIBF)/aeronomars/moldiff.F ; \
	$(AR) r $(LIBO)/libaeronomars.a moldiff.o ; $(RM) moldiff.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(molvis.o) : $(LIBF)/aeronomars/molvis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a molvis.o ; \
	$(COMPILE) $(LIBF)/aeronomars/molvis.F ; \
	$(AR) r $(LIBO)/libaeronomars.a molvis.o ; $(RM) molvis.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(paramfoto_compact.o) : $(LIBF)/aeronomars/paramfoto_compact.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a paramfoto_compact.o ; \
	$(COMPILE) $(LIBF)/aeronomars/paramfoto_compact.F ; \
	$(AR) r $(LIBO)/libaeronomars.a paramfoto_compact.o ; $(RM) paramfoto_compact.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(param_read_e107.o) : $(LIBF)/aeronomars/param_read_e107.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a param_read_e107.o ; \
	$(COMPILE) $(LIBF)/aeronomars/param_read_e107.F ; \
	$(AR) r $(LIBO)/libaeronomars.a param_read_e107.o ; $(RM) param_read_e107.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(param_read.o) : $(LIBF)/aeronomars/param_read.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a param_read.o ; \
	$(COMPILE) $(LIBF)/aeronomars/param_read.F ; \
	$(AR) r $(LIBO)/libaeronomars.a param_read.o ; $(RM) param_read.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(perosat.o) : $(LIBF)/aeronomars/perosat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a perosat.o ; \
	$(COMPILE) $(LIBF)/aeronomars/perosat.F ; \
	$(AR) r $(LIBO)/libaeronomars.a perosat.o ; $(RM) perosat.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(photochemistry.o) : $(LIBF)/aeronomars/photochemistry.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a photochemistry.o ; \
	$(COMPILE90) $(LIBF)/aeronomars/photochemistry.F ; \
	$(AR) r $(LIBO)/libaeronomars.a photochemistry.o ; $(RM) photochemistry.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(read_phototable.o) : $(LIBF)/aeronomars/read_phototable.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/datafile.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a read_phototable.o ; \
	$(COMPILE90) $(LIBF)/aeronomars/read_phototable.F ; \
	$(AR) r $(LIBO)/libaeronomars.a read_phototable.o ; $(RM) read_phototable.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(surfacearea.o) : $(LIBF)/aeronomars/surfacearea.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a surfacearea.o ; \
	$(COMPILE) $(LIBF)/aeronomars/surfacearea.F ; \
	$(AR) r $(LIBO)/libaeronomars.a surfacearea.o ; $(RM) surfacearea.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(thermosphere.o) : $(LIBF)/aeronomars/thermosphere.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v4.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a thermosphere.o ; \
	$(COMPILE) $(LIBF)/aeronomars/thermosphere.F ; \
	$(AR) r $(LIBO)/libaeronomars.a thermosphere.o ; $(RM) thermosphere.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(write_ratecoefficient.o) : $(LIBF)/aeronomars/write_ratecoefficient.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a write_ratecoefficient.o ; \
	$(COMPILE) $(LIBF)/aeronomars/write_ratecoefficient.F ; \
	$(AR) r $(LIBO)/libaeronomars.a write_ratecoefficient.o ; $(RM) write_ratecoefficient.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(calchim.o) : $(LIBF)/aeronomars/calchim.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a calchim.o ; \
	$(COMPILE) $(LIBF)/aeronomars/calchim.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a calchim.o ; $(RM) calchim.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(chemthermos.o) : $(LIBF)/aeronomars/chemthermos.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v4.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a chemthermos.o ; \
	$(COMPILE) $(LIBF)/aeronomars/chemthermos.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a chemthermos.o ; $(RM) chemthermos.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(euvheat.o) : $(LIBF)/aeronomars/euvheat.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v4.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a euvheat.o ; \
	$(COMPILE) $(LIBF)/aeronomars/euvheat.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a euvheat.o ; $(RM) euvheat.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(hcl_deposition.o) : $(LIBF)/aeronomars/hcl_deposition.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a hcl_deposition.o ; \
	$(COMPILE) $(LIBF)/aeronomars/hcl_deposition.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a hcl_deposition.o ; $(RM) hcl_deposition.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(heterogeneous_chlorine.o) : $(LIBF)/aeronomars/heterogeneous_chlorine.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a heterogeneous_chlorine.o ; \
	$(COMPILE) $(LIBF)/aeronomars/heterogeneous_chlorine.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a heterogeneous_chlorine.o ; $(RM) heterogeneous_chlorine.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(inichim_newstart.o) : $(LIBF)/aeronomars/inichim_newstart.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/datafile.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a inichim_newstart.o ; \
	$(COMPTRU90) $(LIBF)/aeronomars/inichim_newstart.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a inichim_newstart.o ; $(RM) inichim_newstart.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(moldiff_red.o) : $(LIBF)/aeronomars/moldiff_red.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h \
$(LIBF)/aeronomars/diffusion.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a moldiff_red.o ; \
	$(COMPILE) $(LIBF)/aeronomars/moldiff_red.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a moldiff_red.o ; $(RM) moldiff_red.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(tgo_profiles.o) : $(LIBF)/aeronomars/tgo_profiles.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a tgo_profiles.o ; \
	$(COMPILE) $(LIBF)/aeronomars/tgo_profiles.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a tgo_profiles.o ; $(RM) tgo_profiles.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(TLM_chem.o) : $(LIBF)/aeronomars/TLM_chem.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a TLM_chem.o ; \
	$(COMPTRU90) $(LIBF)/aeronomars/TLM_chem.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a TLM_chem.o ; $(RM) TLM_chem.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(tlm_clox.o) : $(LIBF)/aeronomars/tlm_clox.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a tlm_clox.o ; \
	$(COMPTRU90) $(LIBF)/aeronomars/tlm_clox.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a tlm_clox.o ; $(RM) tlm_clox.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(tlm_hox.o) : $(LIBF)/aeronomars/tlm_hox.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a tlm_hox.o ; \
	$(COMPTRU90) $(LIBF)/aeronomars/tlm_hox.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a tlm_hox.o ; $(RM) tlm_hox.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(tlm_nox.o) : $(LIBF)/aeronomars/tlm_nox.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a tlm_nox.o ; \
	$(COMPTRU90) $(LIBF)/aeronomars/tlm_nox.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a tlm_nox.o ; $(RM) tlm_nox.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(tlm_ox.o) : $(LIBF)/aeronomars/tlm_ox.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a tlm_ox.o ; \
	$(COMPTRU90) $(LIBF)/aeronomars/tlm_ox.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a tlm_ox.o ; $(RM) tlm_ox.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(tlm_sibem.o) : $(LIBF)/aeronomars/tlm_sibem.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libaeronomars.a tlm_sibem.o ; \
	$(COMPTRU90) $(LIBF)/aeronomars/tlm_sibem.F90 ; \
	$(AR) r $(LIBO)/libaeronomars.a tlm_sibem.o ; $(RM) tlm_sibem.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory bibio
#=======================================================================

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cbrt.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cvmgp.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cvmgt.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(from36.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(gather.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ilsum.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ini36.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(isamax.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(isitor.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ismax.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ismin.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(lnblnk.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxm.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxva.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxv.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mywrite.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(readcoord.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(readstd.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scatter.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scopy.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scopyi.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(sdot.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ssum.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(to36.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(wheneq.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(writestd.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ioipsl_errioipsl.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ioipsl_getincom.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ioipsl_stringop.o)

.PRECIOUS	: $(LIBO)/libbibio.a


# Compilation des membres de la bibliotheque libbibio.a

$(LIBO)/libbibio.a(cbrt.o) : $(LIBF)/bibio/cbrt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a cbrt.o ; \
	$(COMPILE) $(LIBF)/bibio/cbrt.F ; \
	$(AR) r $(LIBO)/libbibio.a cbrt.o ; $(RM) cbrt.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(cvmgp.o) : $(LIBF)/bibio/cvmgp.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a cvmgp.o ; \
	$(COMPILE) $(LIBF)/bibio/cvmgp.F ; \
	$(AR) r $(LIBO)/libbibio.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(cvmgt.o) : $(LIBF)/bibio/cvmgt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a cvmgt.o ; \
	$(COMPILE) $(LIBF)/bibio/cvmgt.F ; \
	$(AR) r $(LIBO)/libbibio.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(from36.o) : $(LIBF)/bibio/from36.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a from36.o ; \
	$(COMPILE) $(LIBF)/bibio/from36.F ; \
	$(AR) r $(LIBO)/libbibio.a from36.o ; $(RM) from36.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(gather.o) : $(LIBF)/bibio/gather.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a gather.o ; \
	$(COMPILE) $(LIBF)/bibio/gather.F ; \
	$(AR) r $(LIBO)/libbibio.a gather.o ; $(RM) gather.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ilsum.o) : $(LIBF)/bibio/ilsum.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ilsum.o ; \
	$(COMPILE) $(LIBF)/bibio/ilsum.F ; \
	$(AR) r $(LIBO)/libbibio.a ilsum.o ; $(RM) ilsum.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ini36.o) : $(LIBF)/bibio/ini36.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ini36.o ; \
	$(COMPILE) $(LIBF)/bibio/ini36.F ; \
	$(AR) r $(LIBO)/libbibio.a ini36.o ; $(RM) ini36.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(isamax.o) : $(LIBF)/bibio/isamax.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a isamax.o ; \
	$(COMPILE) $(LIBF)/bibio/isamax.F ; \
	$(AR) r $(LIBO)/libbibio.a isamax.o ; $(RM) isamax.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(isitor.o) : $(LIBF)/bibio/isitor.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a isitor.o ; \
	$(COMPILE) $(LIBF)/bibio/isitor.F ; \
	$(AR) r $(LIBO)/libbibio.a isitor.o ; $(RM) isitor.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ismax.o) : $(LIBF)/bibio/ismax.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ismax.o ; \
	$(COMPILE) $(LIBF)/bibio/ismax.F ; \
	$(AR) r $(LIBO)/libbibio.a ismax.o ; $(RM) ismax.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ismin.o) : $(LIBF)/bibio/ismin.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ismin.o ; \
	$(COMPILE) $(LIBF)/bibio/ismin.F ; \
	$(AR) r $(LIBO)/libbibio.a ismin.o ; $(RM) ismin.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(lnblnk.o) : $(LIBF)/bibio/lnblnk.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a lnblnk.o ; \
	$(COMPILE) $(LIBF)/bibio/lnblnk.F ; \
	$(AR) r $(LIBO)/libbibio.a lnblnk.o ; $(RM) lnblnk.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxm.o) : $(LIBF)/bibio/mxm.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a mxm.o ; \
	$(COMPILE) $(LIBF)/bibio/mxm.F ; \
	$(AR) r $(LIBO)/libbibio.a mxm.o ; $(RM) mxm.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxva.o) : $(LIBF)/bibio/mxva.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a mxva.o ; \
	$(COMPILE) $(LIBF)/bibio/mxva.F ; \
	$(AR) r $(LIBO)/libbibio.a mxva.o ; $(RM) mxva.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxv.o) : $(LIBF)/bibio/mxv.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a mxv.o ; \
	$(COMPILE) $(LIBF)/bibio/mxv.F ; \
	$(AR) r $(LIBO)/libbibio.a mxv.o ; $(RM) mxv.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mywrite.o) : $(LIBF)/bibio/mywrite.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a mywrite.o ; \
	$(COMPILE) $(LIBF)/bibio/mywrite.F ; \
	$(AR) r $(LIBO)/libbibio.a mywrite.o ; $(RM) mywrite.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(readcoord.o) : $(LIBF)/bibio/readcoord.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a readcoord.o ; \
	$(COMPILE) $(LIBF)/bibio/readcoord.F ; \
	$(AR) r $(LIBO)/libbibio.a readcoord.o ; $(RM) readcoord.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(readstd.o) : $(LIBF)/bibio/readstd.F \
$(LIBF)/bibio/lmdstd.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a readstd.o ; \
	$(COMPILE) $(LIBF)/bibio/readstd.F ; \
	$(AR) r $(LIBO)/libbibio.a readstd.o ; $(RM) readstd.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scatter.o) : $(LIBF)/bibio/scatter.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a scatter.o ; \
	$(COMPILE) $(LIBF)/bibio/scatter.F ; \
	$(AR) r $(LIBO)/libbibio.a scatter.o ; $(RM) scatter.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scopy.o) : $(LIBF)/bibio/scopy.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a scopy.o ; \
	$(COMPILE) $(LIBF)/bibio/scopy.F ; \
	$(AR) r $(LIBO)/libbibio.a scopy.o ; $(RM) scopy.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scopyi.o) : $(LIBF)/bibio/scopyi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a scopyi.o ; \
	$(COMPILE) $(LIBF)/bibio/scopyi.F ; \
	$(AR) r $(LIBO)/libbibio.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(sdot.o) : $(LIBF)/bibio/sdot.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a sdot.o ; \
	$(COMPILE) $(LIBF)/bibio/sdot.F ; \
	$(AR) r $(LIBO)/libbibio.a sdot.o ; $(RM) sdot.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ssum.o) : $(LIBF)/bibio/ssum.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ssum.o ; \
	$(COMPILE) $(LIBF)/bibio/ssum.F ; \
	$(AR) r $(LIBO)/libbibio.a ssum.o ; $(RM) ssum.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(to36.o) : $(LIBF)/bibio/to36.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a to36.o ; \
	$(COMPILE) $(LIBF)/bibio/to36.F ; \
	$(AR) r $(LIBO)/libbibio.a to36.o ; $(RM) to36.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(wheneq.o) : $(LIBF)/bibio/wheneq.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a wheneq.o ; \
	$(COMPILE) $(LIBF)/bibio/wheneq.F ; \
	$(AR) r $(LIBO)/libbibio.a wheneq.o ; $(RM) wheneq.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(writestd.o) : $(LIBF)/bibio/writestd.F \
$(LIBF)/bibio/lmdstd.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a writestd.o ; \
	$(COMPILE) $(LIBF)/bibio/writestd.F ; \
	$(AR) r $(LIBO)/libbibio.a writestd.o ; $(RM) writestd.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ioipsl_errioipsl.o) : $(LIBF)/bibio/ioipsl_errioipsl.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ioipsl_errioipsl.o ; \
	$(COMPTRU90) $(LIBF)/bibio/ioipsl_errioipsl.F90 ; \
	mv $(MOD_LOC_DIR)/ioipsl_errioipsl.$(MOD_SUFFIX) $(LIBO)/ioipsl_errioipsl.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libbibio.a ioipsl_errioipsl.o ; $(RM) ioipsl_errioipsl.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ioipsl_getincom.o) : $(LIBF)/bibio/ioipsl_getincom.F90 \
$(LIBO)/libbibio.a(ioipsl_errioipsl.o) \
$(LIBO)/libbibio.a(ioipsl_stringop.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ioipsl_getincom.o ; \
	$(COMPTRU90) $(LIBF)/bibio/ioipsl_getincom.F90 ; \
	mv $(MOD_LOC_DIR)/ioipsl_getincom.$(MOD_SUFFIX) $(LIBO)/ioipsl_getincom.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libbibio.a ioipsl_getincom.o ; $(RM) ioipsl_getincom.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ioipsl_stringop.o) : $(LIBF)/bibio/ioipsl_stringop.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libbibio.a ioipsl_stringop.o ; \
	$(COMPTRU90) $(LIBF)/bibio/ioipsl_stringop.F90 ; \
	mv $(MOD_LOC_DIR)/ioipsl_stringop.$(MOD_SUFFIX) $(LIBO)/ioipsl_stringop.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libbibio.a ioipsl_stringop.o ; $(RM) ioipsl_stringop.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory dyn3d
#=======================================================================

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(abort_gcm.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(addfi.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(addit.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(advect.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(adv_h2o.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(amont.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(amont_qsat.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(bernoui.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(caldyn0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(caldyn.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(calfis.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(coefpoly.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(convflu.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(convmas.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(coordij.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(covcont.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(covnat.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(defrun_new.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(diagadv.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dissip.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(disvert.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(diverg.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divergf.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(diverg_gam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divergst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divgrad22.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divgrad2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(divgrad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dteta1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dudv1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dudv2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dump2d.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dynetat0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dynredem.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(enercin.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(exner.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(exner_hyb.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(flumass.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fluxstoke.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxhyp.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxy.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxyhyper.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxysinus.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fyhyp.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(geopot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gradiv2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gradiv.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_dyn_fi.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_fi_dyn.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grid_atob.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grid_noro1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(grid_noro.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_int_dyn.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(groupe.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(groupeun.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_u_scal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_v_scal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(heavyside.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniadvtrac.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(ini_archive.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniav.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniconst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(inidissip.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(inigeom.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(inigrads.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniinterp_h.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(initial0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(integrd.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(interp_horiz.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(interp_vert.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien_gam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien_rot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(laplacien_rotgam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(lectflux.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(lect_start_archive.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(massbar.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(massbarxy.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(massdair.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(minmax.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(multipl.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(multscal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgrad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgrad_gam.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgradst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgraro2.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(nxgrarot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(pression.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(psextbar.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(ran1.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(readhead_NC.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotat.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotatf.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotat_nfil.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(rotatst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(scal_wind.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(sortvarc0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(sortvarc.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(spline.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(splint.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(sponge.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(test_period.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(tourpot.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(traceur.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(tracvl.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vanleer.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vitvert.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vlsplt.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vlspltqs.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(wind_scal.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(wrgrads.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(write_archive.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(writediagdyn.o)

.PRECIOUS	: $(LIBO)/libdyn3d.a


# Compilation des membres de la bibliotheque libdyn3d.a

$(LIBO)/libdyn3d.a(abort_gcm.o) : $(LIBF)/dyn3d/abort_gcm.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a abort_gcm.o ; \
	$(COMPILE) $(LIBF)/dyn3d/abort_gcm.F ; \
	$(AR) r $(LIBO)/libdyn3d.a abort_gcm.o ; $(RM) abort_gcm.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(addfi.o) : $(LIBF)/dyn3d/addfi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a addfi.o ; \
	$(COMPILE) $(LIBF)/dyn3d/addfi.F ; \
	$(AR) r $(LIBO)/libdyn3d.a addfi.o ; $(RM) addfi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(addit.o) : $(LIBF)/dyn3d/addit.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a addit.o ; \
	$(COMPILE) $(LIBF)/dyn3d/addit.F ; \
	$(AR) r $(LIBO)/libdyn3d.a addit.o ; $(RM) addit.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(advect.o) : $(LIBF)/dyn3d/advect.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/ener.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a advect.o ; \
	$(COMPILE) $(LIBF)/dyn3d/advect.F ; \
	$(AR) r $(LIBO)/libdyn3d.a advect.o ; $(RM) advect.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(adv_h2o.o) : $(LIBF)/dyn3d/adv_h2o.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a adv_h2o.o ; \
	$(COMPILE) $(LIBF)/dyn3d/adv_h2o.F ; \
	$(AR) r $(LIBO)/libdyn3d.a adv_h2o.o ; $(RM) adv_h2o.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(amont.o) : $(LIBF)/dyn3d/amont.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a amont.o ; \
	$(COMPILE) $(LIBF)/dyn3d/amont.F ; \
	$(AR) r $(LIBO)/libdyn3d.a amont.o ; $(RM) amont.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(amont_qsat.o) : $(LIBF)/dyn3d/amont_qsat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a amont_qsat.o ; \
	$(COMPILE) $(LIBF)/dyn3d/amont_qsat.F ; \
	$(AR) r $(LIBO)/libdyn3d.a amont_qsat.o ; $(RM) amont_qsat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(bernoui.o) : $(LIBF)/dyn3d/bernoui.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a bernoui.o ; \
	$(COMPILE) $(LIBF)/dyn3d/bernoui.F ; \
	$(AR) r $(LIBO)/libdyn3d.a bernoui.o ; $(RM) bernoui.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(caldyn0.o) : $(LIBF)/dyn3d/caldyn0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a caldyn0.o ; \
	$(COMPILE) $(LIBF)/dyn3d/caldyn0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a caldyn0.o ; $(RM) caldyn0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(caldyn.o) : $(LIBF)/dyn3d/caldyn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a caldyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d/caldyn.F ; \
	$(AR) r $(LIBO)/libdyn3d.a caldyn.o ; $(RM) caldyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(calfis.o) : $(LIBF)/dyn3d/calfis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/control.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a calfis.o ; \
	$(COMPILE) $(LIBF)/dyn3d/calfis.F ; \
	$(AR) r $(LIBO)/libdyn3d.a calfis.o ; $(RM) calfis.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(coefpoly.o) : $(LIBF)/dyn3d/coefpoly.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a coefpoly.o ; \
	$(COMPILE) $(LIBF)/dyn3d/coefpoly.F ; \
	$(AR) r $(LIBO)/libdyn3d.a coefpoly.o ; $(RM) coefpoly.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(convflu.o) : $(LIBF)/dyn3d/convflu.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a convflu.o ; \
	$(COMPILE) $(LIBF)/dyn3d/convflu.F ; \
	$(AR) r $(LIBO)/libdyn3d.a convflu.o ; $(RM) convflu.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(convmas.o) : $(LIBF)/dyn3d/convmas.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a convmas.o ; \
	$(COMPILE) $(LIBF)/dyn3d/convmas.F ; \
	$(AR) r $(LIBO)/libdyn3d.a convmas.o ; $(RM) convmas.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(coordij.o) : $(LIBF)/dyn3d/coordij.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a coordij.o ; \
	$(COMPILE) $(LIBF)/dyn3d/coordij.F ; \
	$(AR) r $(LIBO)/libdyn3d.a coordij.o ; $(RM) coordij.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(covcont.o) : $(LIBF)/dyn3d/covcont.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a covcont.o ; \
	$(COMPILE) $(LIBF)/dyn3d/covcont.F ; \
	$(AR) r $(LIBO)/libdyn3d.a covcont.o ; $(RM) covcont.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(covnat.o) : $(LIBF)/dyn3d/covnat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a covnat.o ; \
	$(COMPILE) $(LIBF)/dyn3d/covnat.F ; \
	$(AR) r $(LIBO)/libdyn3d.a covnat.o ; $(RM) covnat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(defrun_new.o) : $(LIBF)/dyn3d/defrun_new.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/comdissnew.h \
$(LIBF)/dyn3d/sponge.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a defrun_new.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/defrun_new.F ; \
	$(AR) r $(LIBO)/libdyn3d.a defrun_new.o ; $(RM) defrun_new.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(diagadv.o) : $(LIBF)/dyn3d/diagadv.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a diagadv.o ; \
	$(COMPILE) $(LIBF)/dyn3d/diagadv.F ; \
	$(AR) r $(LIBO)/libdyn3d.a diagadv.o ; $(RM) diagadv.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dissip.o) : $(LIBF)/dyn3d/dissip.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comdissnew.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dissip.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dissip.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dissip.o ; $(RM) dissip.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(disvert.o) : $(LIBF)/dyn3d/disvert.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a disvert.o ; \
	$(COMPILE) $(LIBF)/dyn3d/disvert.F ; \
	$(AR) r $(LIBO)/libdyn3d.a disvert.o ; $(RM) disvert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(diverg.o) : $(LIBF)/dyn3d/diverg.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a diverg.o ; \
	$(COMPILE) $(LIBF)/dyn3d/diverg.F ; \
	$(AR) r $(LIBO)/libdyn3d.a diverg.o ; $(RM) diverg.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divergf.o) : $(LIBF)/dyn3d/divergf.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divergf.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divergf.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divergf.o ; $(RM) divergf.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(diverg_gam.o) : $(LIBF)/dyn3d/diverg_gam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a diverg_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/diverg_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a diverg_gam.o ; $(RM) diverg_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divergst.o) : $(LIBF)/dyn3d/divergst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divergst.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divergst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divergst.o ; $(RM) divergst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divgrad22.o) : $(LIBF)/dyn3d/divgrad22.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divgrad22.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divgrad22.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divgrad22.o ; $(RM) divgrad22.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divgrad2.o) : $(LIBF)/dyn3d/divgrad2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divgrad2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divgrad2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divgrad2.o ; $(RM) divgrad2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(divgrad.o) : $(LIBF)/dyn3d/divgrad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a divgrad.o ; \
	$(COMPILE) $(LIBF)/dyn3d/divgrad.F ; \
	$(AR) r $(LIBO)/libdyn3d.a divgrad.o ; $(RM) divgrad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dteta1.o) : $(LIBF)/dyn3d/dteta1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dteta1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dteta1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dteta1.o ; $(RM) dteta1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dudv1.o) : $(LIBF)/dyn3d/dudv1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dudv1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dudv1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dudv1.o ; $(RM) dudv1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dudv2.o) : $(LIBF)/dyn3d/dudv2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dudv2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dudv2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dudv2.o ; $(RM) dudv2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dump2d.o) : $(LIBF)/dyn3d/dump2d.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dump2d.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dump2d.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dump2d.o ; $(RM) dump2d.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dynetat0.o) : $(LIBF)/dyn3d/dynetat0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/advtrac.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dynetat0.o ; \
	$(COMPILE90) $(LIBF)/dyn3d/dynetat0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dynetat0.o ; $(RM) dynetat0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dynredem.o) : $(LIBF)/dyn3d/dynredem.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/advtrac.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/advtrac.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a dynredem.o ; \
	$(COMPILE) $(LIBF)/dyn3d/dynredem.F ; \
	$(AR) r $(LIBO)/libdyn3d.a dynredem.o ; $(RM) dynredem.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(enercin.o) : $(LIBF)/dyn3d/enercin.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a enercin.o ; \
	$(COMPILE) $(LIBF)/dyn3d/enercin.F ; \
	$(AR) r $(LIBO)/libdyn3d.a enercin.o ; $(RM) enercin.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(exner.o) : $(LIBF)/dyn3d/exner.F \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a exner.o ; \
	$(COMPILE) $(LIBF)/dyn3d/exner.F ; \
	$(AR) r $(LIBO)/libdyn3d.a exner.o ; $(RM) exner.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(exner_hyb.o) : $(LIBF)/dyn3d/exner_hyb.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a exner_hyb.o ; \
	$(COMPILE) $(LIBF)/dyn3d/exner_hyb.F ; \
	$(AR) r $(LIBO)/libdyn3d.a exner_hyb.o ; $(RM) exner_hyb.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(flumass.o) : $(LIBF)/dyn3d/flumass.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a flumass.o ; \
	$(COMPILE) $(LIBF)/dyn3d/flumass.F ; \
	$(AR) r $(LIBO)/libdyn3d.a flumass.o ; $(RM) flumass.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fluxstoke.o) : $(LIBF)/dyn3d/fluxstoke.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/tracstoke.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fluxstoke.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fluxstoke.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fluxstoke.o ; $(RM) fluxstoke.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxhyp.o) : $(LIBF)/dyn3d/fxhyp.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxhyp.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxhyp.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxhyp.o ; $(RM) fxhyp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxy.o) : $(LIBF)/dyn3d/fxy.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxy_new.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxy.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxy.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxy.o ; $(RM) fxy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxyhyper.o) : $(LIBF)/dyn3d/fxyhyper.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxyhyper.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxyhyper.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxyhyper.o ; $(RM) fxyhyper.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxysinus.o) : $(LIBF)/dyn3d/fxysinus.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxy_sin.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fxysinus.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fxysinus.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fxysinus.o ; $(RM) fxysinus.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fyhyp.o) : $(LIBF)/dyn3d/fyhyp.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a fyhyp.o ; \
	$(COMPILE) $(LIBF)/dyn3d/fyhyp.F ; \
	$(AR) r $(LIBO)/libdyn3d.a fyhyp.o ; $(RM) fyhyp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(geopot.o) : $(LIBF)/dyn3d/geopot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a geopot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/geopot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a geopot.o ; $(RM) geopot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grad.o) : $(LIBF)/dyn3d/grad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grad.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grad.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grad.o ; $(RM) grad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gradiv2.o) : $(LIBF)/dyn3d/gradiv2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gradiv2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gradiv2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gradiv2.o ; $(RM) gradiv2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gradiv.o) : $(LIBF)/dyn3d/gradiv.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gradiv.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gradiv.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gradiv.o ; $(RM) gradiv.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_dyn_fi.o) : $(LIBF)/dyn3d/gr_dyn_fi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_dyn_fi.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_dyn_fi.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_dyn_fi.o ; $(RM) gr_dyn_fi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_fi_dyn.o) : $(LIBF)/dyn3d/gr_fi_dyn.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_fi_dyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_fi_dyn.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_fi_dyn.o ; $(RM) gr_fi_dyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grid_atob.o) : $(LIBF)/dyn3d/grid_atob.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grid_atob.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grid_atob.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grid_atob.o ; $(RM) grid_atob.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grid_noro1.o) : $(LIBF)/dyn3d/grid_noro1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grid_noro1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grid_noro1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grid_noro1.o ; $(RM) grid_noro1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(grid_noro.o) : $(LIBF)/dyn3d/grid_noro.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a grid_noro.o ; \
	$(COMPILE) $(LIBF)/dyn3d/grid_noro.F ; \
	$(AR) r $(LIBO)/libdyn3d.a grid_noro.o ; $(RM) grid_noro.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_int_dyn.o) : $(LIBF)/dyn3d/gr_int_dyn.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_int_dyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_int_dyn.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_int_dyn.o ; $(RM) gr_int_dyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(groupe.o) : $(LIBF)/dyn3d/groupe.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a groupe.o ; \
	$(COMPILE) $(LIBF)/dyn3d/groupe.F ; \
	$(AR) r $(LIBO)/libdyn3d.a groupe.o ; $(RM) groupe.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(groupeun.o) : $(LIBF)/dyn3d/groupeun.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a groupeun.o ; \
	$(COMPILE) $(LIBF)/dyn3d/groupeun.F ; \
	$(AR) r $(LIBO)/libdyn3d.a groupeun.o ; $(RM) groupeun.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_u_scal.o) : $(LIBF)/dyn3d/gr_u_scal.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_u_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_u_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_u_scal.o ; $(RM) gr_u_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_v_scal.o) : $(LIBF)/dyn3d/gr_v_scal.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a gr_v_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/gr_v_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_v_scal.o ; $(RM) gr_v_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(heavyside.o) : $(LIBF)/dyn3d/heavyside.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a heavyside.o ; \
	$(COMPILE) $(LIBF)/dyn3d/heavyside.F ; \
	$(AR) r $(LIBO)/libdyn3d.a heavyside.o ; $(RM) heavyside.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniadvtrac.o) : $(LIBF)/dyn3d/iniadvtrac.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/advtrac.h \
$(LIBF)/dyn3d/control.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a iniadvtrac.o ; \
	$(COMPILE) $(LIBF)/dyn3d/iniadvtrac.F ; \
	$(AR) r $(LIBO)/libdyn3d.a iniadvtrac.o ; $(RM) iniadvtrac.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(ini_archive.o) : $(LIBF)/dyn3d/ini_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/phymars/comsoil.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a ini_archive.o ; \
	$(COMPILE) $(LIBF)/dyn3d/ini_archive.F ; \
	$(AR) r $(LIBO)/libdyn3d.a ini_archive.o ; $(RM) ini_archive.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniav.o) : $(LIBF)/dyn3d/iniav.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comav.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a iniav.o ; \
	$(COMPILE) $(LIBF)/dyn3d/iniav.F ; \
	$(AR) r $(LIBO)/libdyn3d.a iniav.o ; $(RM) iniav.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniconst.o) : $(LIBF)/dyn3d/iniconst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a iniconst.o ; \
	$(COMPILE) $(LIBF)/dyn3d/iniconst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a iniconst.o ; $(RM) iniconst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(inidissip.o) : $(LIBF)/dyn3d/inidissip.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/control.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a inidissip.o ; \
	$(COMPILE) $(LIBF)/dyn3d/inidissip.F ; \
	$(AR) r $(LIBO)/libdyn3d.a inidissip.o ; $(RM) inidissip.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(inigeom.o) : $(LIBF)/dyn3d/inigeom.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comdissnew.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a inigeom.o ; \
	$(COMPILE) $(LIBF)/dyn3d/inigeom.F ; \
	$(AR) r $(LIBO)/libdyn3d.a inigeom.o ; $(RM) inigeom.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(inigrads.o) : $(LIBF)/dyn3d/inigrads.F \
$(LIBF)/dyn3d/gradsdef.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a inigrads.o ; \
	$(COMPILE) $(LIBF)/dyn3d/inigrads.F ; \
	$(AR) r $(LIBO)/libdyn3d.a inigrads.o ; $(RM) inigrads.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniinterp_h.o) : $(LIBF)/dyn3d/iniinterp_h.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a iniinterp_h.o ; \
	$(COMPILE) $(LIBF)/dyn3d/iniinterp_h.F ; \
	$(AR) r $(LIBO)/libdyn3d.a iniinterp_h.o ; $(RM) iniinterp_h.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(initial0.o) : $(LIBF)/dyn3d/initial0.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a initial0.o ; \
	$(COMPILE) $(LIBF)/dyn3d/initial0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a initial0.o ; $(RM) initial0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(integrd.o) : $(LIBF)/dyn3d/integrd.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a integrd.o ; \
	$(COMPILE) $(LIBF)/dyn3d/integrd.F ; \
	$(AR) r $(LIBO)/libdyn3d.a integrd.o ; $(RM) integrd.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(interp_horiz.o) : $(LIBF)/dyn3d/interp_horiz.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a interp_horiz.o ; \
	$(COMPILE) $(LIBF)/dyn3d/interp_horiz.F ; \
	$(AR) r $(LIBO)/libdyn3d.a interp_horiz.o ; $(RM) interp_horiz.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(interp_vert.o) : $(LIBF)/dyn3d/interp_vert.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a interp_vert.o ; \
	$(COMPILE) $(LIBF)/dyn3d/interp_vert.F ; \
	$(AR) r $(LIBO)/libdyn3d.a interp_vert.o ; $(RM) interp_vert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien.o) : $(LIBF)/dyn3d/laplacien.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien.o ; $(RM) laplacien.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien_gam.o) : $(LIBF)/dyn3d/laplacien_gam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien_gam.o ; $(RM) laplacien_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien_rot.o) : $(LIBF)/dyn3d/laplacien_rot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien_rot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien_rot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien_rot.o ; $(RM) laplacien_rot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(laplacien_rotgam.o) : $(LIBF)/dyn3d/laplacien_rotgam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a laplacien_rotgam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/laplacien_rotgam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a laplacien_rotgam.o ; $(RM) laplacien_rotgam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(lectflux.o) : $(LIBF)/dyn3d/lectflux.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/tracstoke.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a lectflux.o ; \
	$(COMPILE) $(LIBF)/dyn3d/lectflux.F ; \
	$(AR) r $(LIBO)/libdyn3d.a lectflux.o ; $(RM) lectflux.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(lect_start_archive.o) : $(LIBF)/dyn3d/lect_start_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comsoil.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/bibio/lmdstd.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/dyn3d/advtrac.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a lect_start_archive.o ; \
	$(COMPILE) $(LIBF)/dyn3d/lect_start_archive.F ; \
	$(AR) r $(LIBO)/libdyn3d.a lect_start_archive.o ; $(RM) lect_start_archive.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(massbar.o) : $(LIBF)/dyn3d/massbar.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a massbar.o ; \
	$(COMPILE) $(LIBF)/dyn3d/massbar.F ; \
	$(AR) r $(LIBO)/libdyn3d.a massbar.o ; $(RM) massbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(massbarxy.o) : $(LIBF)/dyn3d/massbarxy.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a massbarxy.o ; \
	$(COMPILE) $(LIBF)/dyn3d/massbarxy.F ; \
	$(AR) r $(LIBO)/libdyn3d.a massbarxy.o ; $(RM) massbarxy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(massdair.o) : $(LIBF)/dyn3d/massdair.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a massdair.o ; \
	$(COMPILE) $(LIBF)/dyn3d/massdair.F ; \
	$(AR) r $(LIBO)/libdyn3d.a massdair.o ; $(RM) massdair.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(minmax.o) : $(LIBF)/dyn3d/minmax.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a minmax.o ; \
	$(COMPILE) $(LIBF)/dyn3d/minmax.F ; \
	$(AR) r $(LIBO)/libdyn3d.a minmax.o ; $(RM) minmax.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(multipl.o) : $(LIBF)/dyn3d/multipl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a multipl.o ; \
	$(COMPILE) $(LIBF)/dyn3d/multipl.F ; \
	$(AR) r $(LIBO)/libdyn3d.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(multscal.o) : $(LIBF)/dyn3d/multscal.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a multscal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/multscal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a multscal.o ; $(RM) multscal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgrad.o) : $(LIBF)/dyn3d/nxgrad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgrad.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgrad.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgrad.o ; $(RM) nxgrad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgrad_gam.o) : $(LIBF)/dyn3d/nxgrad_gam.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgrad_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgrad_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgrad_gam.o ; $(RM) nxgrad_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgradst.o) : $(LIBF)/dyn3d/nxgradst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgradst.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgradst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgradst.o ; $(RM) nxgradst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgraro2.o) : $(LIBF)/dyn3d/nxgraro2.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgraro2.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgraro2.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgraro2.o ; $(RM) nxgraro2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(nxgrarot.o) : $(LIBF)/dyn3d/nxgrarot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissipn.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a nxgrarot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/nxgrarot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a nxgrarot.o ; $(RM) nxgrarot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(pression.o) : $(LIBF)/dyn3d/pression.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a pression.o ; \
	$(COMPILE) $(LIBF)/dyn3d/pression.F ; \
	$(AR) r $(LIBO)/libdyn3d.a pression.o ; $(RM) pression.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(psextbar.o) : $(LIBF)/dyn3d/psextbar.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a psextbar.o ; \
	$(COMPILE) $(LIBF)/dyn3d/psextbar.F ; \
	$(AR) r $(LIBO)/libdyn3d.a psextbar.o ; $(RM) psextbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(ran1.o) : $(LIBF)/dyn3d/ran1.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a ran1.o ; \
	$(COMPILE) $(LIBF)/dyn3d/ran1.F ; \
	$(AR) r $(LIBO)/libdyn3d.a ran1.o ; $(RM) ran1.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(readhead_NC.o) : $(LIBF)/dyn3d/readhead_NC.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/description.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a readhead_NC.o ; \
	$(COMPILE) $(LIBF)/dyn3d/readhead_NC.F ; \
	$(AR) r $(LIBO)/libdyn3d.a readhead_NC.o ; $(RM) readhead_NC.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotat.o) : $(LIBF)/dyn3d/rotat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotat.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotat.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotat.o ; $(RM) rotat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotatf.o) : $(LIBF)/dyn3d/rotatf.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotatf.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotatf.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotatf.o ; $(RM) rotatf.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotat_nfil.o) : $(LIBF)/dyn3d/rotat_nfil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotat_nfil.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotat_nfil.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotat_nfil.o ; $(RM) rotat_nfil.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(rotatst.o) : $(LIBF)/dyn3d/rotatst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a rotatst.o ; \
	$(COMPILE) $(LIBF)/dyn3d/rotatst.F ; \
	$(AR) r $(LIBO)/libdyn3d.a rotatst.o ; $(RM) rotatst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(scal_wind.o) : $(LIBF)/dyn3d/scal_wind.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a scal_wind.o ; \
	$(COMPILE) $(LIBF)/dyn3d/scal_wind.F ; \
	$(AR) r $(LIBO)/libdyn3d.a scal_wind.o ; $(RM) scal_wind.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(sortvarc0.o) : $(LIBF)/dyn3d/sortvarc0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a sortvarc0.o ; \
	$(COMPILE) $(LIBF)/dyn3d/sortvarc0.F ; \
	$(AR) r $(LIBO)/libdyn3d.a sortvarc0.o ; $(RM) sortvarc0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(sortvarc.o) : $(LIBF)/dyn3d/sortvarc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a sortvarc.o ; \
	$(COMPILE) $(LIBF)/dyn3d/sortvarc.F ; \
	$(AR) r $(LIBO)/libdyn3d.a sortvarc.o ; $(RM) sortvarc.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(spline.o) : $(LIBF)/dyn3d/spline.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a spline.o ; \
	$(COMPILE) $(LIBF)/dyn3d/spline.F ; \
	$(AR) r $(LIBO)/libdyn3d.a spline.o ; $(RM) spline.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(splint.o) : $(LIBF)/dyn3d/splint.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a splint.o ; \
	$(COMPILE) $(LIBF)/dyn3d/splint.F ; \
	$(AR) r $(LIBO)/libdyn3d.a splint.o ; $(RM) splint.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(sponge.o) : $(LIBF)/dyn3d/sponge.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comdissip.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/sponge.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a sponge.o ; \
	$(COMPILE) $(LIBF)/dyn3d/sponge.F ; \
	$(AR) r $(LIBO)/libdyn3d.a sponge.o ; $(RM) sponge.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(test_period.o) : $(LIBF)/dyn3d/test_period.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a test_period.o ; \
	$(COMPILE) $(LIBF)/dyn3d/test_period.F ; \
	$(AR) r $(LIBO)/libdyn3d.a test_period.o ; $(RM) test_period.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(tourpot.o) : $(LIBF)/dyn3d/tourpot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/logic.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a tourpot.o ; \
	$(COMPILE) $(LIBF)/dyn3d/tourpot.F ; \
	$(AR) r $(LIBO)/libdyn3d.a tourpot.o ; $(RM) tourpot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(traceur.o) : $(LIBF)/dyn3d/traceur.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a traceur.o ; \
	$(COMPILE) $(LIBF)/dyn3d/traceur.F ; \
	$(AR) r $(LIBO)/libdyn3d.a traceur.o ; $(RM) traceur.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(tracvl.o) : $(LIBF)/dyn3d/tracvl.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a tracvl.o ; \
	$(COMPILE) $(LIBF)/dyn3d/tracvl.F ; \
	$(AR) r $(LIBO)/libdyn3d.a tracvl.o ; $(RM) tracvl.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vanleer.o) : $(LIBF)/dyn3d/vanleer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vanleer.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vanleer.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vanleer.o ; $(RM) vanleer.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vitvert.o) : $(LIBF)/dyn3d/vitvert.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vitvert.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vitvert.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vitvert.o ; $(RM) vitvert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vlsplt.o) : $(LIBF)/dyn3d/vlsplt.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vlsplt.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vlsplt.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vlsplt.o ; $(RM) vlsplt.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vlspltqs.o) : $(LIBF)/dyn3d/vlspltqs.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a vlspltqs.o ; \
	$(COMPILE) $(LIBF)/dyn3d/vlspltqs.F ; \
	$(AR) r $(LIBO)/libdyn3d.a vlspltqs.o ; $(RM) vlspltqs.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(wind_scal.o) : $(LIBF)/dyn3d/wind_scal.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a wind_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d/wind_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d.a wind_scal.o ; $(RM) wind_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(wrgrads.o) : $(LIBF)/dyn3d/wrgrads.F \
$(LIBF)/dyn3d/gradsdef.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a wrgrads.o ; \
	$(COMPILE) $(LIBF)/dyn3d/wrgrads.F ; \
	$(AR) r $(LIBO)/libdyn3d.a wrgrads.o ; $(RM) wrgrads.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(write_archive.o) : $(LIBF)/dyn3d/write_archive.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/temps.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a write_archive.o ; \
	$(COMPILE) $(LIBF)/dyn3d/write_archive.F ; \
	$(AR) r $(LIBO)/libdyn3d.a write_archive.o ; $(RM) write_archive.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(writediagdyn.o) : $(LIBF)/dyn3d/writediagdyn.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libdyn3d.a writediagdyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d/writediagdyn.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a writediagdyn.o ; $(RM) writediagdyn.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory filtrez
#=======================================================================

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(acc.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(eigen.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(eigen_sort.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(filtreg.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(inifgn.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(inifilr.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(jacobi.o)

.PRECIOUS	: $(LIBO)/libfiltrez.a


# Compilation des membres de la bibliotheque libfiltrez.a

$(LIBO)/libfiltrez.a(acc.o) : $(LIBF)/filtrez/acc.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a acc.o ; \
	$(COMPILE) $(LIBF)/filtrez/acc.F ; \
	$(AR) r $(LIBO)/libfiltrez.a acc.o ; $(RM) acc.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(eigen.o) : $(LIBF)/filtrez/eigen.F \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a eigen.o ; \
	$(COMPILE) $(LIBF)/filtrez/eigen.F ; \
	$(AR) r $(LIBO)/libfiltrez.a eigen.o ; $(RM) eigen.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(eigen_sort.o) : $(LIBF)/filtrez/eigen_sort.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a eigen_sort.o ; \
	$(COMPILE) $(LIBF)/filtrez/eigen_sort.F ; \
	$(AR) r $(LIBO)/libfiltrez.a eigen_sort.o ; $(RM) eigen_sort.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(filtreg.o) : $(LIBF)/filtrez/filtreg.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/filtrez/parafilt.h \
$(LIBF)/filtrez/coefils.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a filtreg.o ; \
	$(COMPILE) $(LIBF)/filtrez/filtreg.F ; \
	$(AR) r $(LIBO)/libfiltrez.a filtreg.o ; $(RM) filtreg.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(inifgn.o) : $(LIBF)/filtrez/inifgn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/filtrez/coefils.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a inifgn.o ; \
	$(COMPILE) $(LIBF)/filtrez/inifgn.F ; \
	$(AR) r $(LIBO)/libfiltrez.a inifgn.o ; $(RM) inifgn.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(inifilr.o) : $(LIBF)/filtrez/inifilr.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/filtrez/parafilt.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/filtrez/coefils.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a inifilr.o ; \
	$(COMPILE) $(LIBF)/filtrez/inifilr.F ; \
	$(AR) r $(LIBO)/libfiltrez.a inifilr.o ; $(RM) inifilr.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(jacobi.o) : $(LIBF)/filtrez/jacobi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libfiltrez.a jacobi.o ; \
	$(COMPILE) $(LIBF)/filtrez/jacobi.F ; \
	$(AR) r $(LIBO)/libfiltrez.a jacobi.o ; $(RM) jacobi.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phymars
#=======================================================================

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(aerave.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(aeropacity.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(aeroptproperties.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(blackl.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(blendrad.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(calldrag_noro.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(callradite.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(callsedim.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(checkfile.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(co2snow.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(convadj.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(cvmgp.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(cvmgt.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(datareadnc.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(drag_noro.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(dustdevil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(dustlift.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(eofdump.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(flusv.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(growthrate.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(gwprofil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(gwstress.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(improvedclouds.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(ineofdump.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(inifis.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(iniorbit.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(inistats.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(initracer.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(iniwrite.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(interp_line.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(ismax.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(ismin.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwb.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwdiff.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwflux.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwmain.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwtt.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwu.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwxb.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwxd.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwxn.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mcd_dust.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mcd_phys.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mcd_trac.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mucorr.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mufract.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(multipl.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(newcondens.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(newsedim.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nirco2abs.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nir_leedat.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nlte_aux.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nlte_calc.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nltecool.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nlte_setup.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nlte_tcool.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nlthermeq.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nuclea.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(orbite.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(orodrag.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(orosetup.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(pbl_parameters.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(phyetat0.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(physdem1.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(physiq.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(profile.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(scatter.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(scopyi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(sig.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(simpleclouds.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(soil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(soil_settings.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(soil_tifeedback.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(solang.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(solarlong.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(sugwd.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(sulw.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(surfini.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(swmain.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(swrayleigh.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(swr_fouquart.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(swr_toon.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(tabfi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(TLM_initialise.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(updatereffrad.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vdif_cd.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vdifc.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vdif_kc.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vlz_fi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(watercloud.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(watersat.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(writediagfi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(writeg1d.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(yamada4.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(zerophys.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(albedocaps.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(calltherm_interface.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(def_var.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(geticecover.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(iniwritesoil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lbfgsb_calls.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mkstat.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(param_slope.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(read_dust_scenario.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(scalar_perturbation.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(suaer.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(thermcell_dqup.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(thermcell_main_mars.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(TLM_ADJ_out.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(TLMvars.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(TLM_vdif.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(updaterad.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(writediagsoil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(wstats.o)

.PRECIOUS	: $(LIBO)/libphymars.a


# Compilation des membres de la bibliotheque libphymars.a

$(LIBO)/libphymars.a(aerave.o) : $(LIBF)/phymars/aerave.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a aerave.o ; \
	$(COMPILE) $(LIBF)/phymars/aerave.F ; \
	$(AR) r $(LIBO)/libphymars.a aerave.o ; $(RM) aerave.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(aeropacity.o) : $(LIBF)/phymars/aeropacity.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/aerkind.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a aeropacity.o ; \
	$(COMPILE90) $(LIBF)/phymars/aeropacity.F ; \
	$(AR) r $(LIBO)/libphymars.a aeropacity.o ; $(RM) aeropacity.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(aeroptproperties.o) : $(LIBF)/phymars/aeroptproperties.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a aeroptproperties.o ; \
	$(COMPILE) $(LIBF)/phymars/aeroptproperties.F ; \
	$(AR) r $(LIBO)/libphymars.a aeroptproperties.o ; $(RM) aeroptproperties.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(blackl.o) : $(LIBF)/phymars/blackl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a blackl.o ; \
	$(COMPILE) $(LIBF)/phymars/blackl.F ; \
	$(AR) r $(LIBO)/libphymars.a blackl.o ; $(RM) blackl.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(blendrad.o) : $(LIBF)/phymars/blendrad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/nlteparams.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a blendrad.o ; \
	$(COMPILE) $(LIBF)/phymars/blendrad.F ; \
	$(AR) r $(LIBO)/libphymars.a blendrad.o ; $(RM) blendrad.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(calldrag_noro.o) : $(LIBF)/phymars/calldrag_noro.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a calldrag_noro.o ; \
	$(COMPILE) $(LIBF)/phymars/calldrag_noro.F ; \
	$(AR) r $(LIBO)/libphymars.a calldrag_noro.o ; $(RM) calldrag_noro.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(callradite.o) : $(LIBF)/phymars/callradite.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/aerkind.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a callradite.o ; \
	$(COMPILE) $(LIBF)/phymars/callradite.F ; \
	$(AR) r $(LIBO)/libphymars.a callradite.o ; $(RM) callradite.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(callsedim.o) : $(LIBF)/phymars/callsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libphymars.a(updaterad.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a callsedim.o ; \
	$(COMPILE90) $(LIBF)/phymars/callsedim.F ; \
	$(AR) r $(LIBO)/libphymars.a callsedim.o ; $(RM) callsedim.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(checkfile.o) : $(LIBF)/phymars/checkfile.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a checkfile.o ; \
	$(COMPILE) $(LIBF)/phymars/checkfile.F ; \
	$(AR) r $(LIBO)/libphymars.a checkfile.o ; $(RM) checkfile.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(co2snow.o) : $(LIBF)/phymars/co2snow.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a co2snow.o ; \
	$(COMPILE) $(LIBF)/phymars/co2snow.F ; \
	$(AR) r $(LIBO)/libphymars.a co2snow.o ; $(RM) co2snow.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(convadj.o) : $(LIBF)/phymars/convadj.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a convadj.o ; \
	$(COMPILE) $(LIBF)/phymars/convadj.F ; \
	$(AR) r $(LIBO)/libphymars.a convadj.o ; $(RM) convadj.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(cvmgp.o) : $(LIBF)/phymars/cvmgp.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a cvmgp.o ; \
	$(COMPILE) $(LIBF)/phymars/cvmgp.F ; \
	$(AR) r $(LIBO)/libphymars.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(cvmgt.o) : $(LIBF)/phymars/cvmgt.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a cvmgt.o ; \
	$(COMPILE) $(LIBF)/phymars/cvmgt.F ; \
	$(AR) r $(LIBO)/libphymars.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(datareadnc.o) : $(LIBF)/phymars/datareadnc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phymars/datafile.h \
$(LIBF)/bibio/lmdstd.h \
$(LIBF)/grid/fxyprim.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a datareadnc.o ; \
	$(COMPILE90) $(LIBF)/phymars/datareadnc.F ; \
	$(AR) r $(LIBO)/libphymars.a datareadnc.o ; $(RM) datareadnc.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(drag_noro.o) : $(LIBF)/phymars/drag_noro.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a drag_noro.o ; \
	$(COMPILE) $(LIBF)/phymars/drag_noro.F ; \
	$(AR) r $(LIBO)/libphymars.a drag_noro.o ; $(RM) drag_noro.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(dustdevil.o) : $(LIBF)/phymars/dustdevil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a dustdevil.o ; \
	$(COMPILE) $(LIBF)/phymars/dustdevil.F ; \
	$(AR) r $(LIBO)/libphymars.a dustdevil.o ; $(RM) dustdevil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(dustlift.o) : $(LIBF)/phymars/dustlift.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a dustlift.o ; \
	$(COMPILE) $(LIBF)/phymars/dustlift.F ; \
	$(AR) r $(LIBO)/libphymars.a dustlift.o ; $(RM) dustlift.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(eofdump.o) : $(LIBF)/phymars/eofdump.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/eofdump.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a eofdump.o ; \
	$(COMPILE) $(LIBF)/phymars/eofdump.F ; \
	$(AR) r $(LIBO)/libphymars.a eofdump.o ; $(RM) eofdump.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(flusv.o) : $(LIBF)/phymars/flusv.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a flusv.o ; \
	$(COMPILE) $(LIBF)/phymars/flusv.F ; \
	$(AR) r $(LIBO)/libphymars.a flusv.o ; $(RM) flusv.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(growthrate.o) : $(LIBF)/phymars/growthrate.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/microphys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a growthrate.o ; \
	$(COMPILE) $(LIBF)/phymars/growthrate.F ; \
	$(AR) r $(LIBO)/libphymars.a growthrate.o ; $(RM) growthrate.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(gwprofil.o) : $(LIBF)/phymars/gwprofil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a gwprofil.o ; \
	$(COMPILE) $(LIBF)/phymars/gwprofil.F ; \
	$(AR) r $(LIBO)/libphymars.a gwprofil.o ; $(RM) gwprofil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(gwstress.o) : $(LIBF)/phymars/gwstress.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a gwstress.o ; \
	$(COMPILE) $(LIBF)/phymars/gwstress.F ; \
	$(AR) r $(LIBO)/libphymars.a gwstress.o ; $(RM) gwstress.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(improvedclouds.o) : $(LIBF)/phymars/improvedclouds.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/microphys.h \
$(LIBF)/aeronomars/conc.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libphymars.a(updaterad.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a improvedclouds.o ; \
	$(COMPILE90) $(LIBF)/phymars/improvedclouds.F ; \
	$(AR) r $(LIBO)/libphymars.a improvedclouds.o ; $(RM) improvedclouds.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(ineofdump.o) : $(LIBF)/phymars/ineofdump.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/eofdump.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a ineofdump.o ; \
	$(COMPILE) $(LIBF)/phymars/ineofdump.F ; \
	$(AR) r $(LIBO)/libphymars.a ineofdump.o ; $(RM) ineofdump.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(inifis.o) : $(LIBF)/phymars/inifis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comsaison.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/datafile.h \
$(LIBF)/phymars/slope.h \
$(LIBF)/phymars/microphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comsoil.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a inifis.o ; \
	$(COMPILE90) $(LIBF)/phymars/inifis.F ; \
	$(AR) r $(LIBO)/libphymars.a inifis.o ; $(RM) inifis.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(iniorbit.o) : $(LIBF)/phymars/iniorbit.F \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a iniorbit.o ; \
	$(COMPILE) $(LIBF)/phymars/iniorbit.F ; \
	$(AR) r $(LIBO)/libphymars.a iniorbit.o ; $(RM) iniorbit.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(inistats.o) : $(LIBF)/phymars/inistats.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phymars/statto.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a inistats.o ; \
	$(COMPILE) $(LIBF)/phymars/inistats.F ; \
	$(AR) r $(LIBO)/libphymars.a inistats.o ; $(RM) inistats.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(initracer.o) : $(LIBF)/phymars/initracer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/dyn3d/advtrac.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a initracer.o ; \
	$(COMPILE) $(LIBF)/phymars/initracer.F ; \
	$(AR) r $(LIBO)/libphymars.a initracer.o ; $(RM) initracer.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(iniwrite.o) : $(LIBF)/phymars/iniwrite.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comsoil.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a iniwrite.o ; \
	$(COMPILE) $(LIBF)/phymars/iniwrite.F ; \
	$(AR) r $(LIBO)/libphymars.a iniwrite.o ; $(RM) iniwrite.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(interp_line.o) : $(LIBF)/phymars/interp_line.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a interp_line.o ; \
	$(COMPILE) $(LIBF)/phymars/interp_line.F ; \
	$(AR) r $(LIBO)/libphymars.a interp_line.o ; $(RM) interp_line.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(ismax.o) : $(LIBF)/phymars/ismax.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a ismax.o ; \
	$(COMPILE) $(LIBF)/phymars/ismax.F ; \
	$(AR) r $(LIBO)/libphymars.a ismax.o ; $(RM) ismax.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(ismin.o) : $(LIBF)/phymars/ismin.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a ismin.o ; \
	$(COMPILE) $(LIBF)/phymars/ismin.F ; \
	$(AR) r $(LIBO)/libphymars.a ismin.o ; $(RM) ismin.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwb.o) : $(LIBF)/phymars/lwb.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwb.o ; \
	$(COMPILE) $(LIBF)/phymars/lwb.F ; \
	$(AR) r $(LIBO)/libphymars.a lwb.o ; $(RM) lwb.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwdiff.o) : $(LIBF)/phymars/lwdiff.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwdiff.o ; \
	$(COMPILE) $(LIBF)/phymars/lwdiff.F ; \
	$(AR) r $(LIBO)/libphymars.a lwdiff.o ; $(RM) lwdiff.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwflux.o) : $(LIBF)/phymars/lwflux.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwflux.o ; \
	$(COMPILE) $(LIBF)/phymars/lwflux.F ; \
	$(AR) r $(LIBO)/libphymars.a lwflux.o ; $(RM) lwflux.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwi.o) : $(LIBF)/phymars/lwi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwi.o ; \
	$(COMPILE) $(LIBF)/phymars/lwi.F ; \
	$(AR) r $(LIBO)/libphymars.a lwi.o ; $(RM) lwi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwmain.o) : $(LIBF)/phymars/lwmain.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwmain.o ; \
	$(COMPILE) $(LIBF)/phymars/lwmain.F ; \
	$(AR) r $(LIBO)/libphymars.a lwmain.o ; $(RM) lwmain.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwtt.o) : $(LIBF)/phymars/lwtt.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwtt.o ; \
	$(COMPILE) $(LIBF)/phymars/lwtt.F ; \
	$(AR) r $(LIBO)/libphymars.a lwtt.o ; $(RM) lwtt.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwu.o) : $(LIBF)/phymars/lwu.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwu.o ; \
	$(COMPILE) $(LIBF)/phymars/lwu.F ; \
	$(AR) r $(LIBO)/libphymars.a lwu.o ; $(RM) lwu.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwxb.o) : $(LIBF)/phymars/lwxb.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwxb.o ; \
	$(COMPILE) $(LIBF)/phymars/lwxb.F ; \
	$(AR) r $(LIBO)/libphymars.a lwxb.o ; $(RM) lwxb.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwxd.o) : $(LIBF)/phymars/lwxd.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwxd.o ; \
	$(COMPILE) $(LIBF)/phymars/lwxd.F ; \
	$(AR) r $(LIBO)/libphymars.a lwxd.o ; $(RM) lwxd.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwxn.o) : $(LIBF)/phymars/lwxn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lwxn.o ; \
	$(COMPILE) $(LIBF)/phymars/lwxn.F ; \
	$(AR) r $(LIBO)/libphymars.a lwxn.o ; $(RM) lwxn.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mcd_dust.o) : $(LIBF)/phymars/mcd_dust.F \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a mcd_dust.o ; \
	$(COMPILE90) $(LIBF)/phymars/mcd_dust.F ; \
	$(AR) r $(LIBO)/libphymars.a mcd_dust.o ; $(RM) mcd_dust.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mcd_phys.o) : $(LIBF)/phymars/mcd_phys.F \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a mcd_phys.o ; \
	$(COMPILE90) $(LIBF)/phymars/mcd_phys.F ; \
	$(AR) r $(LIBO)/libphymars.a mcd_phys.o ; $(RM) mcd_phys.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mcd_trac.o) : $(LIBF)/phymars/mcd_trac.F \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a mcd_trac.o ; \
	$(COMPILE90) $(LIBF)/phymars/mcd_trac.F ; \
	$(AR) r $(LIBO)/libphymars.a mcd_trac.o ; $(RM) mcd_trac.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mucorr.o) : $(LIBF)/phymars/mucorr.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a mucorr.o ; \
	$(COMPILE) $(LIBF)/phymars/mucorr.F ; \
	$(AR) r $(LIBO)/libphymars.a mucorr.o ; $(RM) mucorr.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mufract.o) : $(LIBF)/phymars/mufract.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a mufract.o ; \
	$(COMPILE) $(LIBF)/phymars/mufract.F ; \
	$(AR) r $(LIBO)/libphymars.a mufract.o ; $(RM) mufract.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(multipl.o) : $(LIBF)/phymars/multipl.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a multipl.o ; \
	$(COMPILE) $(LIBF)/phymars/multipl.F ; \
	$(AR) r $(LIBO)/libphymars.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(newcondens.o) : $(LIBF)/phymars/newcondens.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a newcondens.o ; \
	$(COMPILE) $(LIBF)/phymars/newcondens.F ; \
	$(AR) r $(LIBO)/libphymars.a newcondens.o ; $(RM) newcondens.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(newsedim.o) : $(LIBF)/phymars/newsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a newsedim.o ; \
	$(COMPILE) $(LIBF)/phymars/newsedim.F ; \
	$(AR) r $(LIBO)/libphymars.a newsedim.o ; $(RM) newsedim.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nirco2abs.o) : $(LIBF)/phymars/nirco2abs.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/phymars/nirdata.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nirco2abs.o ; \
	$(COMPILE) $(LIBF)/phymars/nirco2abs.F ; \
	$(AR) r $(LIBO)/libphymars.a nirco2abs.o ; $(RM) nirco2abs.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nir_leedat.o) : $(LIBF)/phymars/nir_leedat.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nir_leedat.o ; \
	$(COMPILE) $(LIBF)/phymars/nir_leedat.F ; \
	$(AR) r $(LIBO)/libphymars.a nir_leedat.o ; $(RM) nir_leedat.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nlte_aux.o) : $(LIBF)/phymars/nlte_aux.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nlte_aux.o ; \
	$(COMPILE) $(LIBF)/phymars/nlte_aux.F ; \
	$(AR) r $(LIBO)/libphymars.a nlte_aux.o ; $(RM) nlte_aux.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nlte_calc.o) : $(LIBF)/phymars/nlte_calc.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nlte_calc.o ; \
	$(COMPILE) $(LIBF)/phymars/nlte_calc.F ; \
	$(AR) r $(LIBO)/libphymars.a nlte_calc.o ; $(RM) nlte_calc.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nltecool.o) : $(LIBF)/phymars/nltecool.F \
$(LIBF)/phymars/nltedata.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nltecool.o ; \
	$(COMPILE) $(LIBF)/phymars/nltecool.F ; \
	$(AR) r $(LIBO)/libphymars.a nltecool.o ; $(RM) nltecool.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nlte_setup.o) : $(LIBF)/phymars/nlte_setup.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nlte_setup.o ; \
	$(COMPILE) $(LIBF)/phymars/nlte_setup.F ; \
	$(AR) r $(LIBO)/libphymars.a nlte_setup.o ; $(RM) nlte_setup.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nlte_tcool.o) : $(LIBF)/phymars/nlte_tcool.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nlte_tcool.o ; \
	$(COMPILE) $(LIBF)/phymars/nlte_tcool.F ; \
	$(AR) r $(LIBO)/libphymars.a nlte_tcool.o ; $(RM) nlte_tcool.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nlthermeq.o) : $(LIBF)/phymars/nlthermeq.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/nlteparams.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nlthermeq.o ; \
	$(COMPILE) $(LIBF)/phymars/nlthermeq.F ; \
	$(AR) r $(LIBO)/libphymars.a nlthermeq.o ; $(RM) nlthermeq.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nuclea.o) : $(LIBF)/phymars/nuclea.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/microphys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a nuclea.o ; \
	$(COMPILE) $(LIBF)/phymars/nuclea.F ; \
	$(AR) r $(LIBO)/libphymars.a nuclea.o ; $(RM) nuclea.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(orbite.o) : $(LIBF)/phymars/orbite.F \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a orbite.o ; \
	$(COMPILE) $(LIBF)/phymars/orbite.F ; \
	$(AR) r $(LIBO)/libphymars.a orbite.o ; $(RM) orbite.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(orodrag.o) : $(LIBF)/phymars/orodrag.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a orodrag.o ; \
	$(COMPILE) $(LIBF)/phymars/orodrag.F ; \
	$(AR) r $(LIBO)/libphymars.a orodrag.o ; $(RM) orodrag.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(orosetup.o) : $(LIBF)/phymars/orosetup.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a orosetup.o ; \
	$(COMPILE) $(LIBF)/phymars/orosetup.F ; \
	$(AR) r $(LIBO)/libphymars.a orosetup.o ; $(RM) orosetup.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(pbl_parameters.o) : $(LIBF)/phymars/pbl_parameters.F \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a pbl_parameters.o ; \
	$(COMPILE) $(LIBF)/phymars/pbl_parameters.F ; \
	$(AR) r $(LIBO)/libphymars.a pbl_parameters.o ; $(RM) pbl_parameters.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(phyetat0.o) : $(LIBF)/phymars/phyetat0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/dyn3d/advtrac.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a phyetat0.o ; \
	$(COMPILE90) $(LIBF)/phymars/phyetat0.F ; \
	$(AR) r $(LIBO)/libphymars.a phyetat0.o ; $(RM) phyetat0.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(physdem1.o) : $(LIBF)/phymars/physdem1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comdissnew.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/advtrac.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/clesph0.h \
$(LIBF)/grid/fxyprim.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comsoil.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a physdem1.o ; \
	$(COMPILE) $(LIBF)/phymars/physdem1.F ; \
	$(AR) r $(LIBO)/libphymars.a physdem1.o ; $(RM) physdem1.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(physiq.o) : $(LIBF)/phymars/physiq.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comsoil.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comsaison.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/nlteparams.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v4.h \
$(LIBF)/aeronomars/conc.h \
$(LIBF)/phymars/slope.h \
$(LIBF)/dyn3d/advtrac.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a physiq.o ; \
	$(COMPILE90) $(LIBF)/phymars/physiq.F ; \
	$(AR) r $(LIBO)/libphymars.a physiq.o ; $(RM) physiq.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(profile.o) : $(LIBF)/phymars/profile.F \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a profile.o ; \
	$(COMPILE90) $(LIBF)/phymars/profile.F ; \
	$(AR) r $(LIBO)/libphymars.a profile.o ; $(RM) profile.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(scatter.o) : $(LIBF)/phymars/scatter.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a scatter.o ; \
	$(COMPILE) $(LIBF)/phymars/scatter.F ; \
	$(AR) r $(LIBO)/libphymars.a scatter.o ; $(RM) scatter.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(scopyi.o) : $(LIBF)/phymars/scopyi.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a scopyi.o ; \
	$(COMPILE) $(LIBF)/phymars/scopyi.F ; \
	$(AR) r $(LIBO)/libphymars.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(sig.o) : $(LIBF)/phymars/sig.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a sig.o ; \
	$(COMPILE) $(LIBF)/phymars/sig.F ; \
	$(AR) r $(LIBO)/libphymars.a sig.o ; $(RM) sig.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(simpleclouds.o) : $(LIBF)/phymars/simpleclouds.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBO)/libphymars.a(updaterad.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a simpleclouds.o ; \
	$(COMPILE90) $(LIBF)/phymars/simpleclouds.F ; \
	$(AR) r $(LIBO)/libphymars.a simpleclouds.o ; $(RM) simpleclouds.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(soil.o) : $(LIBF)/phymars/soil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comsoil.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a soil.o ; \
	$(COMPILE) $(LIBF)/phymars/soil.F ; \
	$(AR) r $(LIBO)/libphymars.a soil.o ; $(RM) soil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(soil_settings.o) : $(LIBF)/phymars/soil_settings.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comsoil.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a soil_settings.o ; \
	$(COMPILE90) $(LIBF)/phymars/soil_settings.F ; \
	$(AR) r $(LIBO)/libphymars.a soil_settings.o ; $(RM) soil_settings.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(soil_tifeedback.o) : $(LIBF)/phymars/soil_tifeedback.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comsoil.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a soil_tifeedback.o ; \
	$(COMPILE) $(LIBF)/phymars/soil_tifeedback.F ; \
	$(AR) r $(LIBO)/libphymars.a soil_tifeedback.o ; $(RM) soil_tifeedback.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(solang.o) : $(LIBF)/phymars/solang.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a solang.o ; \
	$(COMPILE) $(LIBF)/phymars/solang.F ; \
	$(AR) r $(LIBO)/libphymars.a solang.o ; $(RM) solang.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(solarlong.o) : $(LIBF)/phymars/solarlong.F \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a solarlong.o ; \
	$(COMPILE) $(LIBF)/phymars/solarlong.F ; \
	$(AR) r $(LIBO)/libphymars.a solarlong.o ; $(RM) solarlong.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(sugwd.o) : $(LIBF)/phymars/sugwd.F \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a sugwd.o ; \
	$(COMPILE) $(LIBF)/phymars/sugwd.F ; \
	$(AR) r $(LIBO)/libphymars.a sugwd.o ; $(RM) sugwd.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(sulw.o) : $(LIBF)/phymars/sulw.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a sulw.o ; \
	$(COMPILE) $(LIBF)/phymars/sulw.F ; \
	$(AR) r $(LIBO)/libphymars.a sulw.o ; $(RM) sulw.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(surfini.o) : $(LIBF)/phymars/surfini.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/datafile.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a surfini.o ; \
	$(COMPILE90) $(LIBF)/phymars/surfini.F ; \
	$(AR) r $(LIBO)/libphymars.a surfini.o ; $(RM) surfini.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(swmain.o) : $(LIBF)/phymars/swmain.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a swmain.o ; \
	$(COMPILE) $(LIBF)/phymars/swmain.F ; \
	$(AR) r $(LIBO)/libphymars.a swmain.o ; $(RM) swmain.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(swrayleigh.o) : $(LIBF)/phymars/swrayleigh.F \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a swrayleigh.o ; \
	$(COMPILE) $(LIBF)/phymars/swrayleigh.F ; \
	$(AR) r $(LIBO)/libphymars.a swrayleigh.o ; $(RM) swrayleigh.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(swr_fouquart.o) : $(LIBF)/phymars/swr_fouquart.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a swr_fouquart.o ; \
	$(COMPILE) $(LIBF)/phymars/swr_fouquart.F ; \
	$(AR) r $(LIBO)/libphymars.a swr_fouquart.o ; $(RM) swr_fouquart.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(swr_toon.o) : $(LIBF)/phymars/swr_toon.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a swr_toon.o ; \
	$(COMPILE) $(LIBF)/phymars/swr_toon.F ; \
	$(AR) r $(LIBO)/libphymars.a swr_toon.o ; $(RM) swr_toon.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(tabfi.o) : $(LIBF)/phymars/tabfi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comsoil.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a tabfi.o ; \
	$(COMPILE) $(LIBF)/phymars/tabfi.F ; \
	$(AR) r $(LIBO)/libphymars.a tabfi.o ; $(RM) tabfi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(TLM_initialise.o) : $(LIBF)/phymars/TLM_initialise.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a TLM_initialise.o ; \
	$(COMPILE90) $(LIBF)/phymars/TLM_initialise.F ; \
	$(AR) r $(LIBO)/libphymars.a TLM_initialise.o ; $(RM) TLM_initialise.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(updatereffrad.o) : $(LIBF)/phymars/updatereffrad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/aerkind.h \
$(LIBF)/phymars/yomaer.h \
$(LIBO)/libphymars.a(updaterad.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a updatereffrad.o ; \
	$(COMPILE90) $(LIBF)/phymars/updatereffrad.F ; \
	$(AR) r $(LIBO)/libphymars.a updatereffrad.o ; $(RM) updatereffrad.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vdif_cd.o) : $(LIBF)/phymars/vdif_cd.F \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a vdif_cd.o ; \
	$(COMPILE) $(LIBF)/phymars/vdif_cd.F ; \
	$(AR) r $(LIBO)/libphymars.a vdif_cd.o ; $(RM) vdif_cd.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vdifc.o) : $(LIBF)/phymars/vdifc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/microphys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a vdifc.o ; \
	$(COMPILE90) $(LIBF)/phymars/vdifc.F ; \
	$(AR) r $(LIBO)/libphymars.a vdifc.o ; $(RM) vdifc.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vdif_kc.o) : $(LIBF)/phymars/vdif_kc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a vdif_kc.o ; \
	$(COMPILE) $(LIBF)/phymars/vdif_kc.F ; \
	$(AR) r $(LIBO)/libphymars.a vdif_kc.o ; $(RM) vdif_kc.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vlz_fi.o) : $(LIBF)/phymars/vlz_fi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a vlz_fi.o ; \
	$(COMPILE) $(LIBF)/phymars/vlz_fi.F ; \
	$(AR) r $(LIBO)/libphymars.a vlz_fi.o ; $(RM) vlz_fi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(watercloud.o) : $(LIBF)/phymars/watercloud.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o) \
$(LIBO)/libphymars.a(updaterad.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a watercloud.o ; \
	$(COMPILE90) $(LIBF)/phymars/watercloud.F ; \
	$(AR) r $(LIBO)/libphymars.a watercloud.o ; $(RM) watercloud.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(watersat.o) : $(LIBF)/phymars/watersat.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a watersat.o ; \
	$(COMPILE) $(LIBF)/phymars/watersat.F ; \
	$(AR) r $(LIBO)/libphymars.a watersat.o ; $(RM) watersat.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(writediagfi.o) : $(LIBF)/phymars/writediagfi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a writediagfi.o ; \
	$(COMPILE90) $(LIBF)/phymars/writediagfi.F ; \
	$(AR) r $(LIBO)/libphymars.a writediagfi.o ; $(RM) writediagfi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(writeg1d.o) : $(LIBF)/phymars/writeg1d.F \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a writeg1d.o ; \
	$(COMPILE) $(LIBF)/phymars/writeg1d.F ; \
	$(AR) r $(LIBO)/libphymars.a writeg1d.o ; $(RM) writeg1d.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(yamada4.o) : $(LIBF)/phymars/yamada4.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a yamada4.o ; \
	$(COMPILE) $(LIBF)/phymars/yamada4.F ; \
	$(AR) r $(LIBO)/libphymars.a yamada4.o ; $(RM) yamada4.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(zerophys.o) : $(LIBF)/phymars/zerophys.F
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a zerophys.o ; \
	$(COMPILE) $(LIBF)/phymars/zerophys.F ; \
	$(AR) r $(LIBO)/libphymars.a zerophys.o ; $(RM) zerophys.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(albedocaps.o) : $(LIBF)/phymars/albedocaps.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/datafile.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a albedocaps.o ; \
	$(COMPTRU90) $(LIBF)/phymars/albedocaps.F90 ; \
	$(AR) r $(LIBO)/libphymars.a albedocaps.o ; $(RM) albedocaps.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(calltherm_interface.o) : $(LIBF)/phymars/calltherm_interface.F90 \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a calltherm_interface.o ; \
	$(COMPTRU90) $(LIBF)/phymars/calltherm_interface.F90 ; \
	$(AR) r $(LIBO)/libphymars.a calltherm_interface.o ; $(RM) calltherm_interface.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(def_var.o) : $(LIBF)/phymars/def_var.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a def_var.o ; \
	$(COMPILE) $(LIBF)/phymars/def_var.F90 ; \
	$(AR) r $(LIBO)/libphymars.a def_var.o ; $(RM) def_var.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(geticecover.o) : $(LIBF)/phymars/geticecover.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a geticecover.o ; \
	$(COMPILE) $(LIBF)/phymars/geticecover.F90 ; \
	$(AR) r $(LIBO)/libphymars.a geticecover.o ; $(RM) geticecover.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(iniwritesoil.o) : $(LIBF)/phymars/iniwritesoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/phymars/comsoil.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a iniwritesoil.o ; \
	$(COMPILE) $(LIBF)/phymars/iniwritesoil.F90 ; \
	$(AR) r $(LIBO)/libphymars.a iniwritesoil.o ; $(RM) iniwritesoil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lbfgsb_calls.o) : $(LIBF)/phymars/lbfgsb_calls.F90 \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/tracer.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a lbfgsb_calls.o ; \
	$(COMPTRU90) $(LIBF)/phymars/lbfgsb_calls.F90 ; \
	$(AR) r $(LIBO)/libphymars.a lbfgsb_calls.o ; $(RM) lbfgsb_calls.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mkstat.o) : $(LIBF)/phymars/mkstat.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/statto.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a mkstat.o ; \
	$(COMPILE) $(LIBF)/phymars/mkstat.F90 ; \
	$(AR) r $(LIBO)/libphymars.a mkstat.o ; $(RM) mkstat.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(param_slope.o) : $(LIBF)/phymars/param_slope.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a param_slope.o ; \
	$(COMPILE) $(LIBF)/phymars/param_slope.F90 ; \
	$(AR) r $(LIBO)/libphymars.a param_slope.o ; $(RM) param_slope.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(read_dust_scenario.o) : $(LIBF)/phymars/read_dust_scenario.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/datafile.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a read_dust_scenario.o ; \
	$(COMPTRU90) $(LIBF)/phymars/read_dust_scenario.F90 ; \
	$(AR) r $(LIBO)/libphymars.a read_dust_scenario.o ; $(RM) read_dust_scenario.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(scalar_perturbation.o) : $(LIBF)/phymars/scalar_perturbation.F90 \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a scalar_perturbation.o ; \
	$(COMPTRU90) $(LIBF)/phymars/scalar_perturbation.F90 ; \
	$(AR) r $(LIBO)/libphymars.a scalar_perturbation.o ; $(RM) scalar_perturbation.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(suaer.o) : $(LIBF)/phymars/suaer.F90 \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/datafile.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/aerkind.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a suaer.o ; \
	$(COMPILE) $(LIBF)/phymars/suaer.F90 ; \
	$(AR) r $(LIBO)/libphymars.a suaer.o ; $(RM) suaer.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(thermcell_dqup.o) : $(LIBF)/phymars/thermcell_dqup.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a thermcell_dqup.o ; \
	$(COMPILE) $(LIBF)/phymars/thermcell_dqup.F90 ; \
	$(AR) r $(LIBO)/libphymars.a thermcell_dqup.o ; $(RM) thermcell_dqup.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(thermcell_main_mars.o) : $(LIBF)/phymars/thermcell_main_mars.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a thermcell_main_mars.o ; \
	$(COMPILE) $(LIBF)/phymars/thermcell_main_mars.F90 ; \
	$(AR) r $(LIBO)/libphymars.a thermcell_main_mars.o ; $(RM) thermcell_main_mars.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(TLM_ADJ_out.o) : $(LIBF)/phymars/TLM_ADJ_out.F90 \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h \
$(LIBO)/libbibio.a(ioipsl_getincom.o)
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a TLM_ADJ_out.o ; \
	$(COMPTRU90) $(LIBF)/phymars/TLM_ADJ_out.F90 ; \
	$(AR) r $(LIBO)/libphymars.a TLM_ADJ_out.o ; $(RM) TLM_ADJ_out.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(TLMvars.o) : $(LIBF)/phymars/TLMvars.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a TLMvars.o ; \
	$(COMPTRU90) $(LIBF)/phymars/TLMvars.F90 ; \
	mv $(MOD_LOC_DIR)/TLMvars.$(MOD_SUFFIX) $(LIBO)/TLMvars.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphymars.a TLMvars.o ; $(RM) TLMvars.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(TLM_vdif.o) : $(LIBF)/phymars/TLM_vdif.F90
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a TLM_vdif.o ; \
	$(COMPTRU90) $(LIBF)/phymars/TLM_vdif.F90 ; \
	$(AR) r $(LIBO)/libphymars.a TLM_vdif.o ; $(RM) TLM_vdif.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(updaterad.o) : $(LIBF)/phymars/updaterad.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a updaterad.o ; \
	$(COMPTRU90) $(LIBF)/phymars/updaterad.F90 ; \
	mv $(MOD_LOC_DIR)/updaterad.$(MOD_SUFFIX) $(LIBO)/updaterad.$(MOD_SUFFIX) ; \
	$(AR) r $(LIBO)/libphymars.a updaterad.o ; $(RM) updaterad.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(writediagsoil.o) : $(LIBF)/phymars/writediagsoil.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/phymars/comsoil.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a writediagsoil.o ; \
	$(COMPILE) $(LIBF)/phymars/writediagsoil.F90 ; \
	$(AR) r $(LIBO)/libphymars.a writediagsoil.o ; $(RM) writediagsoil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(wstats.o) : $(LIBF)/phymars/wstats.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/statto.h
	cd $(LOCAL_DIR); \
	$(AR) d $(LIBO)/libphymars.a wstats.o ; \
	$(COMPILE) $(LIBF)/phymars/wstats.F90 ; \
	$(AR) r $(LIBO)/libphymars.a wstats.o ; $(RM) wstats.o ; \
	cd $(GCM)


