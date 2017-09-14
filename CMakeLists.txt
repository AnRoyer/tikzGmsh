cmake_minimum_required(VERSION 3.5)

project(tikz)

file(GLOB FILES *.cpp)

include_directories(include/gmsh)

link_directories(lib)
set(EXECUTABLE_OUTPUT_PATH ../bin)


add_executable(tikz ${FILES})


target_link_libraries(tikz Gmsh lapack blas)
