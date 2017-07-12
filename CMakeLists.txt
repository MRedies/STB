cmake_minimum_required(VERSION 2.6)
project(STB)

file(GLOB_RECURSE sources  src/*.f90 
                           src/*.F90
                           src/*.h 
                           thirdparty/config_fortran/m_config.f90
                           thirdparty/NPY-for-Fortran/src/endian_swap.f90
                           thirdparty/triangulation/tri_wrapper.f90
                           thirdparty/NPY-for-Fortran/src/npy.F90)



link_directories($ENV{LAPACK_LIB} /bgsys/local/lib thirdparty/triangulation/)
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
    set(dialect "-ffree-form -Wall -fimplicit-none -m64  -m64 -I$ENV{MKLROOT}/include  -L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
    set(bounds "-fbounds-check -Wall")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-std08 -warn all -g -qopenmp -implicitnone -mkl -heap_arrays")
    #set(dialect "-std08 -warn all -g -implicitnone -mkl:sequential -heap_arrays")
    set(bounds "-check all")
    set(debug "-traceback -debug ") 
endif() 
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
	set(bounds "-C")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "IBM")
    MESSAGE(STATUS "IBM")
    target_link_libraries(stb.x lapack esslbg)
    set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
    #set(CMAKE_Fortran_COMPILER scorep)
    #set(dialect "mpixlf2008_r -O3 -qstrict -qarch=qp -qtune=qp -WF,-DIBM_COMPILER_USED -qnullterm")
    set(CMAKE_Fortran_COMPILER mpixlf2008_r)
    set(dialect "-O3 -qstrict -qarch=qp -qtune=qp -WF,-DIBM_COMPILER_USED -qnullterm")
    set(debug "-C")
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
