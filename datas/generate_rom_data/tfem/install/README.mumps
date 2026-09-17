The MUMPS package can be obtained from http://graal.ens-lyon.fr/MUMPS/.

The MUMPS package contains a sparse direct solver and can be compiled in
two ways:
  - sequential version, which is a free alternative for MA41 and MA57.
  - parallel version using MPI.

The mumps add-on contains an interface to the MUMPS package and has been
tested and verified to be working with the following installation of parallel
MUMPS and its dependencies:
  - MUMPS 5.1.2
  - Open MPI 3.0.1
  - Metis 5.1.0
  - Scotch 6.0.5a
  - MKL 18.0.0

Metis and Scotch can be ommited but are strongly recommended in order to
increase performance (time and memory issues). Additional gains can be
expected when OpenMP is also used. A sample Makefile.inc for MUMPS can
be found in the tfem folder install (in the root of tfem). 

You should modify the LIBBLAS macro in the Makefile.inc to either:
  - empty, in which case the LIBBLAS is defined by Mdefs.mk,
  - identical to the LIBBLAS in Mdefs.mk.
otherwise you might link to a system library blas, which is usually slower.
In de first case, the MUMPS library will be built, but the examples will not.


Installation
============

- Adjust the root Mdefs.mk of tfem to use the compiler that is used to
  compile MUMPS.

- Edit Mdefs.mk in the root of tfem by adding the contents of Mdefs_mumps.mk to 
  the end of Mdefs.mk file.

- Adjust Mdefs.mk according to the directions given in the comments.

- In case gfortran is used, remove the first line (ifeq) and the third
  line (endif) statement in the file
         addons/io_utils/Mdefs_io_utils.mk
  since MUMPS will modify the FC variable.

- build tfem using 

    make all mumps

  for a usual build + mumps.

