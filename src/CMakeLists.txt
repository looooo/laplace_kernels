# build library
add_library(laplaceKern2D SHARED laplaceKern2D.cpp ${headers})
add_library(laplaceKern3D SHARED laplaceKern3D.cpp ${headers})

# installing
install(TARGETS laplaceKern2D laplaceKern3D LIBRARY DESTINATION lib)
install(FILES laplaceKern2D.h laplaceKern3D.h DESTINATION include/laplace)

# testing
if (MAKE_TEST)
    add_library(vtkwriter SHARED vtkWriter.cpp)
    target_link_libraries(vtkwriter ${VTK_LIBRARIES})

    add_executable(test2D test2D.cpp)
    target_link_libraries(test2D laplaceKern2D vtkwriter)

    add_executable(test3D test3D.cpp)
    target_link_libraries(test3D laplaceKern3D vtkwriter)
endif()

if (MAKE_PYTHON_BINDINGS)
    add_library(laplaceKern SHARED laplaceKern_pyext.cpp)
    set_target_properties(laplaceKern PROPERTIES PREFIX "")
    target_link_libraries(laplaceKern laplaceKern2D laplaceKern3D)
	set_target_properties(laplaceKern PROPERTIES INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
    install(TARGETS laplaceKern DESTINATION ${python_libs})
endif()