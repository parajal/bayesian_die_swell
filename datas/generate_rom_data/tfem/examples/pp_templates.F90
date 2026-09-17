#if SKIT2
#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on sk_solve', &
    ' - install and compile sparskit2', &
    ' - set preprocessing macro SKIT2 in Mdefs.mk'
end
#endif

#if METIS5
#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on metis5', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro METIS5 in Mdefs.mk'
end
#endif

#if SKIT2 && METIS5
#else
  print '(4(a/),a)', &
    'To run this example:', &
    ' - compile add-ons sk_solve and metis5', &
    ' - install and compile sparskit2', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro SKIT2 and METIS5 in Mdefs.mk'
end
#endif

#if PARDISO
#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on pardiso', &
    ' - add the MKL library for linking with pardiso', &
    ' - set preprocessing macro PARDISO in Mdefs.mk'
end
#endif

#if HSL_EXTRA
#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on hsl_extra', &
    ' - add the libhsl3 library for linking', &
    ' - set preprocessing macro HSL_EXTRA in Mdefs.mk'
end
#endif

#if PARTRACK
#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on particle_tracking', &
    ' - add the nr library for linking', &
    ' - set preprocessing macro PARTRACK in Mdefs.mk'
end
#endif

#if MKL_GMRES && METIS5
#else
  print '(4(a/),a)', &
    'To run this example:', &
    ' - compile add-ons mkl_gmres and metis5', &
    ' - add the MKL library for linking with mkl_gmres', &
    ' - add metis5 lib for linking', &
    ' - set preprocessing macro MKL_GMRES and METIS5 in Mdefs.mk'
end
#endif

#if UMFPACK
#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on umfpack', &
    ' - add the UMFPACK library for linking', &
    ' - set preprocessing macro UMFPACK in Mdefs.mk'
end
#endif

#if FEAST
#else
  print '(3(a/),a)', &
    'To run this example:', &
    ' - compile add-on feast', &
    ' - add the feast library for linking', &
    ' - set preprocessing macro FEAST in Mdefs.mk'
end
#endif
