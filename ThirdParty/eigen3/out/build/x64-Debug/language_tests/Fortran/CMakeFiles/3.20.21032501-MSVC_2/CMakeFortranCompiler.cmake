set(CMAKE_Fortran_COMPILER "C:/Strawberry/c/bin/gfortran.exe")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "4.9.2")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "MinGW")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "C:/Strawberry/c/bin/ar.exe")
set(CMAKE_Fortran_COMPILER_AR "C:/Strawberry/c/bin/gcc-ar.exe")
set(CMAKE_RANLIB "C:/Strawberry/c/bin/ranlib.exe")
set(CMAKE_Fortran_COMPILER_RANLIB "C:/Strawberry/c/bin/gcc-ranlib.exe")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW 1)
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "4")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "C:/Strawberry/c/lib/gcc/i686-w64-mingw32/4.9.2/include;C:/Strawberry/c/include;C:/Strawberry/c/lib/gcc/i686-w64-mingw32/4.9.2/include-fixed;C:/Strawberry/c/i686-w64-mingw32/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;mingw32;gcc_s;gcc;moldname;mingwex;quadmath;m;mingw32;gcc_s;gcc;moldname;mingwex;pthread;advapi32;shell32;user32;kernel32;iconv;mingw32;gcc_s;gcc;moldname;mingwex")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "C:/Strawberry/c/lib/gcc/i686-w64-mingw32/4.9.2;C:/Strawberry/c/lib/gcc;C:/Strawberry/c/i686-w64-mingw32/lib;C:/Strawberry/c/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
