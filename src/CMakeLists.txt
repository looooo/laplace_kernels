# build library
add_library(laplaceKern2D SHARED laplaceKern2D.cpp ${headers})
add_library(laplaceKern3D SHARED laplaceKern3D.cpp ${headers})
add_library(vtkwriter SHARED vtkWriter.cpp)
target_link_libraries(vtkwriter ${VTK_LIBRARIES})

# installing
install(TARGETS laplaceKern2D laplaceKern3D vtkwriter LIBRARY DESTINATION lib/laplace)
install(FILES laplaceKern2D.h laplaceKern3D.h vtkWriter.h DESTINATION include/laplace)

# testing
MESSAGE(${MAKE_TEST})
if (${MAKE_TEST})
    add_executable(test2D test2D.cpp vtkWriter.cpp)
    target_link_libraries(test2D laplaceKern2D vtkwriter)

    add_executable(test3D test3D.cpp vtkWriter.cpp)
    target_link_libraries(test3D laplaceKern3D vtkwriter)
endif()