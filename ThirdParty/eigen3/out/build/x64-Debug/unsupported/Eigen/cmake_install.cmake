# Install script for directory: E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/AdolcForward"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/AlignedVector3"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/ArpackSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/AutoDiff"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/BVH"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/EulerAngles"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/FFT"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/IterativeSolvers"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/KroneckerProduct"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/LevenbergMarquardt"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/MatrixFunctions"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/MoreVectorization"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/MPRealSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/NonLinearOptimization"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/NumericalDiff"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/OpenGLSupport"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/Polynomials"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/Skyline"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/SparseExtra"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/SpecialFunctions"
    "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("E:/Cowain_Algorithm_library/CoVisionLibs/ThirdParty/eigen3/out/build/x64-Debug/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

