set(SRC_FILES
      src/modules/constants.f90       ## libSUFR
      src/modules/date_and_time.f90   ## libSUFR
      src/modules/dummy_variables.f90 ## libSUFR
      src/modules/date_and_time.f90   ## libSUFR
      src/modules/getopt.f90          ## libSUFR
      src/modules/kinds.f90           ## libSUFR
      src/modules/system.f90          ## libSUFR
      src/modules/text.f90            ## libSUFR
      src/modules/commonmod.f90
      src/modules/conductance.f90
      src/modules/gen_multpol_mat.f90
      src/modules/memory_use.f90
      src/modules/mo_kind.f90
      src/modules/mo_quicksort.f90
      src/modules/mo_utils.f90
      src/modules/readmod.f90
      src/modules/writemod.f90
      src/compute/basis2atom_map.f90
      src/compute/c2ci.f90
      src/compute/d2dmo.f90
      src/compute/dip2dipmo.f90
      src/compute/EDOs.f90
      src/compute/FIOs.f90
      src/FIOs-EDOs-mkl.f90
      src/matlib/diasym.f90
      src/matlib/matinv.f90
      src/matlib/mtxm_mkl.f90
      src/matlib/mxm_mkl.f90
      src/matlib/mxmt_mkl.f90
      src/matlib/mtxmxm_mkl.f90
      src/matlib/mxmxmt_mkl.f90
      src/matlib/mxvxmt_mkl.f90
      src/read/read_pol.f90
      src/read/readfchk.f90
      src/read/readfchk2.f90
      src/read/readg09out.f90
      src/sort/indexxabs.f90
   )
