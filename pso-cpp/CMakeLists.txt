cmake_minimum_required(VERSION 3.5)
project(pso)

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(${SRC_DIR})

add_executable(run_experiments ${SRC_DIR}/experiment_runner.cpp ${SOURCES})
