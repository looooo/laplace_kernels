cmake_minimum_required(VERSION 2.6)
project(laplacekernels)

add_executable(laplacekernels main.cpp)

install(TARGETS laplacekernels RUNTIME DESTINATION bin)
