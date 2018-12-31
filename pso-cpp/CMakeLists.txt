cmake_minimum_required(VERSION 3.5)
project(pso)

set(SRC_DIR "${PROJECT_SOURCE_DIR}"/src)

set(CMAKE_CXX_STANDARD 17)

include_directories(${SRC_DIR})

set(SOURCES
  ${SRC_DIR}/particle.cpp
  ${SRC_DIR}/simulation.cpp
  ${SRC_DIR}/experiment.cpp
)

add_executable(run_experiments ${SRC_DIR}/experiment_runner.cpp ${SOURCES})
