#
# Edit for running MUMPS:
#
# 1) Adjust the MUMPSPATH macro to represent your MUMPS installation
# 2) Adjust the LSCOTCHDIR and METISPATH macros.

MUMPSPATH = 

ifdef MUMPSPATH

# MUMPS

  include $(MUMPSPATH)/Makefile.inc # NOTE: redefines macro FC

  LPORDDIR = $(MUMPSPATH)/PORD/lib/
  LPORD    = -L$(LPORDDIR) -lpord

  LMETIS    = -L$(METISPATH) -lmetis

  SCOTCHDIR  = ${HOME}/scotch_6.0.5a
  ISCOTCH    = -I$(SCOTCHDIR)/include
  LSCOTCH    = -L$(SCOTCHDIR)/lib -lesmumps -lscotch -lscotcherr

  LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)

  ifdef LIBSEQNEEDED

    # sequential

    INC = -I$(MUMPSPATH)/libseq -I$(MUMPSPATH)/include
    LIB = -L$(MUMPSPATH)/lib -ldmumps -lmumps_common -L$(MUMPSPATH)/libseq -lmpiseq

  else

    # parallel

    INC = $(INCPAR) -I$(MUMPSPATH)/include $(ISCOTCH)
    LIB = -L$(MUMPSPATH)/lib -ldmumps -lmumps_common $(LIBPAR)

  endif

  FMODINC = -I$(MODPATHINC) $(FMODINCP) $(INC)

  LIBS = $(LIBSP) -L$(TFEMPATH)/lib -ltfem $(LIB) $(LIBBLAS) $(LIBOTHERS) \
         $(LORDERINGS)

endif
