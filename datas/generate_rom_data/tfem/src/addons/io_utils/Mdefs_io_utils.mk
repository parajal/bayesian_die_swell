ifeq ($(FC),gfortran) 
  FFLAGS = $(FCHECK) $(FOPT) $(FDEBUG) $(FMISC) -std=gnu
endif
