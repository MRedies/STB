cmake_minimum_required(VERSION 2.6)
project(STB)

file(GLOB_RECURSE sources  src/*.f90 
                           src/*.h 
                           thirdparty/config_fortran/m_config.f90
                           thirdparty/NPY-for-Fortran/src/npy.F90)
add_executable(stb.x ${sources})

enable_language(Fortran)
set(CMAKE_Fortran_COMPILER mpifort)
set(CMAKE_Fortran_COMPILER_ID "Intel")
#set(CMAKE_Fortran_COMPILER gfortran)
#set(CMAKE_Fortran_COMPILER_ID "GNU")

MESSAGE( STATUS "MKLROOT: $ENV{MKLROOT}")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -Wall -fimplicit-none -m64 -I$ENV{MKLROOT}/include")
    set (CMAKE_EXE_LINKER_FLAGS "-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
    set(bounds "-fbounds-check -Wall")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand=f08 -warn all -implicitnone -g -mkl:sequential -heap_arrays")
    add_definitions(-DINTEL_COMPILER_USED) 
    set(bounds "-check all")
    set(debug "-traceback -debug ") 
endif() 
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
	set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug}")
MESSAGE( STATUS "DEBUG Flag ${CMAKE_Fortran_FLAGS_DEBUG}")


set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
MESSAGE( STATUS "Flags:  ${CMAKE_Fortran_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE( STATUS "cmake_module_path:    " ${CMAKE_MODULE_PATH})

#find_package(MKL REQUIRED)
#if(MKL_FOUND)
    #MESSAGE( STATUS "MKL FOUND")
    #include_directories(${MKL_INCLUDE_DIRS})
    #target_link_libraries (stb.x ${MKL_MINIMAL_LIBRARY})
#else()
    #MESSAGE( STATUS "MKL NOT FOUND")
#endif()
MESSAGE(STATUS "MKL_THREADING_LIBNAME " ${MKL_THREADING_LIBNAME})
#
# Compile.
#
