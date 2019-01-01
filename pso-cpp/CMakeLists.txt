cmake_minimum_required(VERSION 3.9)
project(pso)

find_package(OpenMP REQUIRED)

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 ")

include_directories(${SRC_DIR})

add_executable(run_experiments ${SRC_DIR}/experiment_runner.cpp ${SOURCES})
target_link_libraries(run_experiments PRIVATE OpenMP::OpenMP_CXX)
