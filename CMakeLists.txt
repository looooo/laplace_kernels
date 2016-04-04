cmake_minimum_required(VERSION 2.6)
project(laplacekernels)


ADD_DEFINITIONS(-std=c++11)

find_package(Eigen3)
include_directories( ${EIGEN3_INCLUDE_DIR})

if (${MAKE_TEST})
	find_package(VTK REQUIRED)
	find_package(Qt5Widgets REQUIRED)
	include(${VTK_USE_FILE})
endif()

add_subdirectory(src)
