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


# Open MPI support
option(ENABLE_MPI    "Enable MPI support?" ON)
if(ENABLE_MPI)
	find_package(MPI)
	if(MPI_FOUND)
		list(APPEND COMMON_DEFINITIONS -DHAVE_MPI)
		list(APPEND COMMON_C_FLAGS ${MPI_C_FLAGS})
		list(APPEND COMMON_LIBRARIES ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
		list(APPEND COMMON_INCLUDE_PATHS ${MPI_C_INCLUDE_PATH})
	endif()
endif()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS} ${COMMON_DEFINITIONS}")
add_executable(ps4 ps4.c util.c)
target_link_libraries(ps4 ${COMMON_LIBRARIES})

