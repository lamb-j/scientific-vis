cmake_minimum_required(VERSION 2.8)

PROJECT(isosurface)
#SET(VTK_DIR /usr/local/Cellar/vtk/8.1.0/lib/cmake/vtk-8.1) 

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(isosurface main.cxx isosurface_prefix.cxx isosurface_parallel.cxx tricase.cxx render_surface.cxx)

set(CMAKE_CXX_FLAGS "-std=c++11 -I/usr/local/cuda-8.0/include/")
SET(CMAKE_CXX_LINK_FLAGS "-lOpenCL")

if(VTK_LIBRARIES)
  target_link_libraries(isosurface ${VTK_LIBRARIES})
else()
  target_link_libraries(isosurface vtkHybrid)
endif()
