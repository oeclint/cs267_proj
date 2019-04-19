# - Try to find Eigen2
# Once done, this will define
#
#  Eigen2_FOUND - system has Eigen2
#  Eigen2_INCLUDE_DIR - the Eigen2 include directory

# Include dir
find_path(Eigen3_INCLUDE_DIR
  NAMES eigen3/Eigen/Core
  PATH_SUFFIXES eigen3
)

set(Eigen3_FOUND "NO")
if(Eigen3_INCLUDE_DIR)
  set(Eigen3_FOUND "YES")
endif(Eigen3_INCLUDE_DIR)

