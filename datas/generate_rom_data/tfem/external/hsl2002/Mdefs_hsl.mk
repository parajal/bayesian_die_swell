MODPATHWRITE=.
LIBCREATE=$(TFEMPATH)/lib/lib$(FSOLVERLIB).a
ifeq ($(FC),nagfor) 
  FFLAGS = $(FOPT) $(FMISC) $(FCHECK) -dusty
else
  FFLAGS = $(FOPT) $(FMISC)
endif
