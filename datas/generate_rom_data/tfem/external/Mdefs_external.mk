LIBCREATE=$(TFEMPATH)/lib/lib$(FSOLVERLIB).a
ifeq ($(FC),nagfor) 
  FFLAGS = $(FOPT) $(FMISC) $(FCHECK) -dusty
else ifeq ($(FC),gfortran)
  FFLAGS = $(FOPT) $(FMISC) $(FCHECK) -std=legacy
else
  FFLAGS = $(FOPT) $(FMISC)
endif
