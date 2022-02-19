cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

string(APPEND CMAKE_C_FLAGS_DEBUG " -DDEBUG")
execute_process (
  COMMAND git describe --abbrev=4 --dirty --always --tags
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_definitions(-DVERSION="${GIT_COMMIT_HASH}")
add_definitions(-DN_DIMS=3)

# set the project name
project(symd)

# empty ID = Apple Clang for some reason
if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
set(APPEND CMAKE_C_FLAGS_RELEASE "-O2 -g")
else()
set(CMAKE_C_FLAGS_RELEASE "-O3 -fno-tree-pta -g -march=native -funroll-loops -fmodulo-sched -fmodulo-sched-allow-regmoves")
#set(CMAKE_C_FLAGS_RELEASE "-O3 -g -Wall")
endif()



message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")
message("Compiler is ${CMAKE_CXX_COMPILER_ID}")


# add sources
include_directories(include)
set(SOURCE_FILES src/cJSON.c
                src/harmonic_force.c
                src/lj_force.c
                src/soft_force.c
                src/main_loop.c
                src/box.c
                src/nlist.c
                src/vverlet_integrate.c
                src/util.c
                src/anderson_thermostat.c
                src/bussi_thermostat.c
                src/group.c
                )

# add the executable
add_executable(symd ${SOURCE_FILES})

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(symd GSL::gsl GSL::gslcblas)
target_link_libraries(symd m)
#find_package(OpenMP REQUIRED)
#target_link_libraries(symd OpenMP::OpenMP_CXX)