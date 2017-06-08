cmake_minimum_required(VERSION 2.6)
project(STB)

file(GLOB_RECURSE sources  src/*.f90 
                           src/*.F90
                           src/*.h 
                           thirdparty/config_fortran/m_config.f90
                           thirdparty/NPY-for-Fortran/src/npy.F90)

add_executable(stb.x ${sources})

enable_language(Fortran)
#set(CMAKE_Fortran_COMPILER mpiifort)
#set(CMAKE_Fortran_COMPILER_ID "Intel")
#set(CMAKE_Fortran_COMPILER mpifort)
#set(CMAKE_Fortran_COMPILER_ID "GNU")
#set(CMAKE_Fortran_COMPILER mpifort)
set(CMAKE_Fortran_COMPILER mpixlf2003)
set(CMAKE_Fortran_COMPILER_ID "IBM")

find_package (LAPACK REQUIRED)
target_link_libraries (stb.x ${LAPACK_LIBRARIES})

MESSAGE( STATUS "MKLROOT: $ENV{MKLROOT}")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -Wall -fimplicit-none -m64")
    set(bounds "-fbounds-check -Wall")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand=f08 -warn all -implicitnone -g -heap_arrays")
    set(bounds "-check all")
    set(debug "-traceback -debug ") 
endif() 
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
	set(bounds "-C")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "IBM")
    MESSAGE(STATUS "IBM")
    set(dialect "-O2 -qarch=qp -qtune=qp -ufmt_littleendian={23}")
    set(CMAKE_EXE_LINKER_FLAGS "-qmkshrobj")
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-Wl,-export-dynamic")
    set(debug "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug}")
MESSAGE( STATUS "DEBUG Flag ${CMAKE_Fortran_FLAGS_DEBUG}")


set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
MESSAGE( STATUS "Flags:  ${CMAKE_Fortran_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE( STATUS "cmake_module_path:    " ${CMAKE_MODULE_PATH})

#
# Compile.
#
