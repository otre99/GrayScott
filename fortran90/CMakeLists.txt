cmake_minimum_required (VERSION 2.8)
project (gray_scott_fortran)
enable_language (Fortran)

SET(EXP "user_code" CACHE STRING "User specified experiment")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fopenmp") 
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

SET(SRC src/gs_types.f90 
        src/gs_utils.f90
        src/gs_solver.f90
        src/gs_io.f90
        src/gs_main.f90
        )

add_executable(${PROJECT_NAME} ${SRC})
install (TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .)
