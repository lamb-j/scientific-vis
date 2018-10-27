cmake_minimum_required(VERSION 3.10)

PROJECT(astro)
SET(VTK_DIR /usr/local/Cellar/vtk/8.1.0/lib/cmake/vtk-8.1) 

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(astro astro)

set(CMAKE_CXX_FLAGS "-Ofast")
SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(astro ${VTK_LIBRARIES})
else()
  target_link_libraries(astro vtkHybrid)
endif()
