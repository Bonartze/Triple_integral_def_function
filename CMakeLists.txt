cmake_minimum_required(VERSION 3.25)
project(untitled51 C)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

add_executable(untitled51 main.c integral.h integral.c)
