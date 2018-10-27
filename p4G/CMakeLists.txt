cmake_minimum_required(VERSION 3.10)

PROJECT(proj4G)
SET(VTK_DIR /usr/local/Cellar/vtk/8.1.0/lib/cmake/vtk-8.1) 

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj4G proj4G)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj4G ${VTK_LIBRARIES})
else()
  target_link_libraries(proj4G vtkHybrid)
endif()
