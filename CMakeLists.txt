cmake_minimum_required(VERSION 2.6)
project(laplacekernels)

ADD_DEFINITIONS(-std=c++14)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
SET(CMAKE_BUILD_TYPE REALEASE)
SET(CMAKE_CXX_FLAGS_REALEASE "-O3")

find_package(Eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})

option(MAKE_TEST "test require vtk" ON)
option(MAKE_PYTHON_BINDINGS "create python bindings" ON)

if (MAKE_TEST)
    find_package(VTK REQUIRED)
    find_package(Qt5Widgets REQUIRED)
    include(${VTK_USE_FILE})

    message("create test executable")
endif()


if (MAKE_PYTHON_BINDINGS)
    find_package(PythonInterp)
    find_package(PythonLibs)
    EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c 
    "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE python_libs OUTPUT_STRIP_TRAILING_WHITESPACE )
    include_directories( ${PYTHON_INCLUDE_DIRS})
endif()

add_subdirectory(src)
