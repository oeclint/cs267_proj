set(CMAKE_CXX_COMPILER "/opt/cray/pe/craype/2.5.15/bin/CC")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "Intel")
set(CMAKE_CXX_COMPILER_VERSION "18.0.1.20171018")
set(CMAKE_CXX_COMPILER_WRAPPER "CrayPrgEnv")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "98")
set(CMAKE_CXX_COMPILE_FEATURES "")
set(CMAKE_CXX98_COMPILE_FEATURES "")
set(CMAKE_CXX11_COMPILE_FEATURES "")
set(CMAKE_CXX14_COMPILE_FEATURES "")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "")
set(CMAKE_CXX_SIMULATE_VERSION "")

set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_LINKER "/usr/common/software/altd/2.0/bin/ld")
set(CMAKE_COMPILER_IS_GNUCXX )
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;mm;CPP)
set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()




set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "fmpich;mpichcxx;darshan;darshan-stubs;z;AtpSigHandler;AtpSigHCommData;pthread;mpichcxx_intel;rt;ugni;pthread;pmi;imf;m;dl;sci_intel_mpi;sci_intel;imf;m;dl;mpich_intel;rt;ugni;pthread;pmi;imf;m;dl;pmi;pthread;alpslli;pthread;wlm_detect;alpsutil;pthread;rca;xpmem;ugni;pthread;udreg;sci_intel;imf;m;pthread;dl;hugetlbfs;stdc++;imf;m;ifcore;ifport;pthread;imf;svml;irng;stdc++;m;ipgo;decimal;stdc++;irc;svml;c;irc_s;dl;c")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/opt/cray/pe/libsci/18.07.1/INTEL/16.0/x86_64/lib;/opt/cray/dmapp/default/lib64;/opt/cray/pe/mpt/7.7.3/gni/mpich-intel/16.0/lib;/usr/common/software/darshan/3.1.4/lib;/opt/cray/rca/2.2.18-6.0.7.1_5.47__g2aa4f39.ari/lib64;/opt/cray/alps/6.6.43-6.0.7.1_5.45__ga796da32.ari/lib64;/opt/cray/xpmem/2.2.15-6.0.7.1_5.11__g7549d06.ari/lib64;/opt/cray/pe/pmi/5.0.14/lib64;/opt/cray/ugni/6.0.14.0-6.0.7.1_3.13__gea11d3d.ari/lib64;/opt/cray/udreg/2.3.2-6.0.7.1_5.13__g5196236.ari/lib64;/opt/cray/pe/atp/2.1.3/libApp;/lib64;/opt/cray/wlm_detect/1.3.3-6.0.7.1_5.6__g7109084.ari/lib64;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64;/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64;/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin;/usr/lib64/gcc/x86_64-suse-linux/4.8;/usr/lib64;/usr/x86_64-suse-linux/lib;/lib;/usr/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
