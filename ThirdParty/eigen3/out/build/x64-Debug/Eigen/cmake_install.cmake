# Install script for directory: E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/out/install/x64-Debug")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "E:/python-pycharm-opencv/MinGW/bin/objdump.exe")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Cholesky"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/CholmodSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Core"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Dense"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Eigen"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Eigenvalues"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Geometry"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Householder"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/IterativeLinearSolvers"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Jacobi"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/LU"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/MetisSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/OrderingMethods"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/PaStiXSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/PardisoSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/QR"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/QtAlignedMalloc"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/SPQRSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/SVD"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/Sparse"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/SparseCholesky"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/SparseCore"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/SparseLU"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/SparseQR"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/StdDeque"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/StdList"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/StdVector"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/SuperLUSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

