cmake_minimum_required(VERSION 3.20.0)
project(fold LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
add_compile_options(-w)
file(GLOB SOURCES src/protein.cpp src/*.cu)
add_executable(fold-bin src/main.cpp ${SOURCES})
include_directories(${EIGEN3_INCLUDE_DIR} src/fold.cuh)

set(PYTHON_EXECUTABLE /usr/bin/python)
find_package(pybind11 CONFIG REQUIRED)
# find_package(PythonLibs REQUIRED)
pybind11_add_module(fold src/pybind.cpp ${SOURCES})
target_include_directories(fold PRIVATE ${pybind11_INCLUDE_DIRS})
