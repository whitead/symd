cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS_RELEASE "-Ofast -march=native -fmodulo-sched-allow-regmoves")

# set the project name
project(simple-md)

# add sources
include_directories(include)
set(SOURCE_FILES src/cJSON.c
                src/harmonic_force.c
                src/lj_force.c
                src/soft_force.c
                src/main_loop.c
                src/min_image.c
                src/nlist.c
                src/vverlet_integrate.c
                src/util.c
                src/anderson_thermostat.c
                src/bussi_thermostat.c
                src/group.c
                )

# add the executable
add_executable(simple-md ${SOURCE_FILES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(simple-md GSL::gsl GSL::gslcblas)
target_link_libraries(simple-md m)
find_package(OpenMP REQUIRED)
target_link_libraries(simple-md OpenMP::OpenMP_CXX)