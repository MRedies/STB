cmake_minimum_required(VERSION 3.5)
project(STB)

set(MKL_MULTI_THREADED FALSE)

enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
	#set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand=f08 -implicitnone")
	#set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
	#set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE( STATUS "cmake_module_path:    " ${CMAKE_MODULE_PATH})

file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(stb.x ${sources})

find_package(MKL REQUIRED)
if(MKL_FOUND)
    MESSAGE( STATUS "MKL FOUND")
    include_directories(${MKL_INCLUDE_DIRS})
    target_link_libraries (stb.x ${MKL_MINIMAL_LIBRARY})
else()
    MESSAGE( STATUS "MKL NOT FOUND")
endif()
#
# Compile.
#
