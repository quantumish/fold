cmake_minimum_required(VERSION 3.10.0)

project(fold)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

if (PYTHON)
    find_package(pybind11 CONFIG REQUIRED)
    include_directories(${pybind11_INCLUDE_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DPYTHON")
    pybind11_add_module(fold fold.cpp)
endif()

if (CPP)
    add_executable(fold ./fold.cpp)
endif()
  
