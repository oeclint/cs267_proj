set(CMAKE_C_COMPILER "/opt/cray/pe/craype/2.5.15/bin/cc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "Intel")
set(CMAKE_C_COMPILER_VERSION "18.0.1.20171018")
set(CMAKE_C_COMPILER_WRAPPER "CrayPrgEnv")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "90")
set(CMAKE_C_COMPILE_FEATURES "")
set(CMAKE_C90_COMPILE_FEATURES "")
set(CMAKE_C99_COMPILE_FEATURES "")
set(CMAKE_C11_COMPILE_FEATURES "")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_SIMULATE_VERSION "")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/common/software/altd/2.0/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()




set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "fmpich;mpichcxx;darshan;darshan-stubs;z;AtpSigHandler;AtpSigHCommData;pthread;sci_intel_mpi;sci_intel;imf;m;dl;mpich_intel;rt;ugni;pthread;pmi;imf;m;dl;pmi;pthread;alpslli;pthread;wlm_detect;alpsutil;pthread;rca;xpmem;ugni;pthread;udreg;sci_intel;imf;m;pthread;dl;hugetlbfs;imf;m;ifcore;ifport;pthread;imf;svml;irng;m;ipgo;decimal;irc;svml;c;irc_s;dl;c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/cray/pe/libsci/18.07.1/INTEL/16.0/x86_64/lib;/opt/cray/dmapp/default/lib64;/opt/cray/pe/mpt/7.7.3/gni/mpich-intel/16.0/lib;/usr/common/software/darshan/3.1.4/lib;/opt/cray/rca/2.2.18-6.0.7.1_5.47__g2aa4f39.ari/lib64;/opt/cray/alps/6.6.43-6.0.7.1_5.45__ga796da32.ari/lib64;/opt/cray/xpmem/2.2.15-6.0.7.1_5.11__g7549d06.ari/lib64;/opt/cray/pe/pmi/5.0.14/lib64;/opt/cray/ugni/6.0.14.0-6.0.7.1_3.13__gea11d3d.ari/lib64;/opt/cray/udreg/2.3.2-6.0.7.1_5.13__g5196236.ari/lib64;/opt/cray/pe/atp/2.1.3/libApp;/lib64;/opt/cray/wlm_detect/1.3.3-6.0.7.1_5.6__g7109084.ari/lib64;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin;/usr/lib64/gcc/x86_64-suse-linux/4.8;/usr/lib64;/usr/x86_64-suse-linux/lib;/lib;/usr/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
