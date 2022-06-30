cmake_minimum_required(VERSION 2.6)
project(STB)

SET(STB_USE_64BIT OFF CACHE BOOL "Use 64-bit interface")

file(GLOB_RECURSE sources  src/*.f90
                           src/*.F90
                           src/*.h
                           thirdparty/config_fortran/m_config.f90
                           thirdparty/triangulation/tri_wrapper.f90
                           thirdparty/NPY-for-Fortran/src/npy.f90)

link_directories($ENV{LAPACK_LIB} /bgsys/local/lib thirdparty/triangulation/ /opt/ibmcmp/xlsmp/bg/3.1/bglib64 $ENV{MKLROOT}/lib/intel64)
enable_language(Fortran)
add_executable(stb.x ${sources})

if(PLATFORM MATCHES "INTEL")
    set(CMAKE_Fortran_COMPILER_ID "Intel")
    set(CMAKE_Fortran_COMPILER mpiifort)
endif()
if(PLATFORM MATCHES "GNU")
    set(CMAKE_Fortran_COMPILER_ID "GNU")
    set(CMAKE_Fortran_COMPILER mpifort)
endif()
if(PLATFORM MATCHES "IBM")
    set(CMAKE_Fortran_COMPILER_ID "IBM")
endif()

MESSAGE(STATUS "Using ${CMAKE_Fortran_COMPILER_ID} compiler")

target_link_libraries (stb.x ${LAPACK_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/triangulation/triangle.o)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(stb.x PRIVATE -C
        -Wall
        -Wextra
        -Wconversion
        -Wline-truncation
        -Wcharacter-truncation
        -Wsurprising
        -Waliasing
        -Wunused-parameter
        -Wno-maybe-uninitialized
        -fwhole-file
        -fimplicit-none
        -fallow-argument-mismatch
        -fbacktrace
        -ffree-form
        -ffree-line-length-0
        -fcheck=all
        -finit-real=nan
        -ftrapv
        -cpp
        -g
        -O0
    )
    set(dialect "-ffree-form -Wall -fimplicit-none -m64  -m64 -I$ENV{MKLROOT}/include  -L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
    set(bounds "-fbounds-check -Wall")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(stb.x PRIVATE -C
        -CB
        -check uninit
        -debug all
        -warn all
        -e18
        -ftrapuv
        -traceback
        -qopenmp
        -fpp
        -g
        -O0
    )
    #set(dialect "-std08 -warn all -g -init=Nan -implicitnone -heap_arrays")
    set(dialect "-std08 -warn all -g -implicitnone -heap_arrays")
    set(bounds "-check all")
    set(debug "-traceback -debug ")
    if(STB_USE_64BIT)
        target_link_libraries(stb.x mkl_intel_ilp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
        set(dialect "${dialect} -i8 -I$ENV{MKLROOT}/include")
        MESSAGE(STATUS "using 64-bit interface")
    else()
        target_link_libraries(stb.x mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
        set(dialect "${dialect} -I$ENV{MKLROOT}/include")
    endif()
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug}")
MESSAGE( STATUS "DEBUG Flag ${CMAKE_Fortran_FLAGS_DEBUG}")


set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
MESSAGE( STATUS "Flags:  ${CMAKE_Fortran_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE( STATUS "cmake_module_path:    " ${CMAKE_MODULE_PATH})

MESSAGE(STATUS ${CMAKE_EXE_LINKER_FLAGS})
MESSAGE(STATUS ${CMAKE_MODULE_LINKER_FLAGS})
MESSAGE(STATUS ${CMAKE_SHARED_LINKER_FLAGS})
MESSAGE(STATUS ${CMAKE_STATIC_LINKER_FLAGS})
#
# Compile.
#
