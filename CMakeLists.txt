cmake_minimum_required(VERSION 3.10)

# set the project name
project(simple-md)

# add sources
include_directories(include)
set(SOURCE_FILES src/cJSON.c
                src/harmonic_force.c
                src/lj_force.c
                src/main_loop.c
                src/min_image.c
                src/nlist.c
                src/vverlet_integrate.c
                src/util.c
                )

# add the executable
add_executable(simple-md ${SOURCE_FILES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(simple-md GSL::gsl GSL::gslcblas)
target_link_libraries(simple-md m)
# find_package(OpenMP REQUIRED)
# target_link_libraries(simple-md OpenMP::OpenMP_CXX)