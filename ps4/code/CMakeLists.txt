cmake_minimum_required(VERSION 2.6)
project(ps4)
enable_language(C)
# The version number
set(ps4_VERSION_MAJOR 0)
set(ps4_VERSION_MINOR 1)



# Open MP support
option(ENABLE_OPENMP "Enable OpenMP support?" ON)
if(ENABLE_OPENMP)
	find_package(OpenMP)
  	if(OPENMP_FOUND)
    		list(APPEND COMMON_DEFINITIONS -DHAVE_OPENMP)
		list(APPEND COMMON_C_FLAGS ${OpenMP_C_FLAGS})
	endif()
endif(ENABLE_OPENMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS} ${COMMON_DEFINITIONS}")
add_executable(ps4 ps4.c util.c)

