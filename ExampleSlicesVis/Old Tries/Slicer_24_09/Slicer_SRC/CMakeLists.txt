cmake_minimum_required(VERSION 2.8)
 
PROJECT(ReadImageData)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(ReadImageData MACOSX_BUNDLE ReadImageData)
 
if(VTK_LIBRARIES)
  target_link_libraries(ReadImageData myVtkInteractorStyleImage ${VTK_LIBRARIES})
else()
  target_link_libraries(ReadImageData myVtkInteractorStyleImage  vtkHybrid vtkWidgets)
endif()