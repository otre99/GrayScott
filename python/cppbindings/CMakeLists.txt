cmake_minimum_required(VERSION 3.14)

project(gray_scott_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -O3")

find_package(Python 3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# set output name
exec_program(${Python_EXECUTABLE}-config ARGS --extension-suffix OUTPUT_VARIABLE LIBSUFFIX)
set(LIBNAME gray_scott_cpp${LIBSUFFIX})
get_filename_component(LIBNAME ${LIBNAME}  NAME_WLE)
message(STATUS ${LIBNAME})

set(SRC main.cpp gs_solver.cpp gs_solver.h)
add_library(${LIBNAME} MODULE ${SRC})
target_link_libraries(${LIBNAME} PUBLIC pybind11::module pthread tbb)

## Or embed the Python interpreter into an executable
#add_executable(myexe main.cpp)
#target_link_libraries(myexe PUBLIC pybind11::embed)
