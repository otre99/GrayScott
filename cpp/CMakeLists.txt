cmake_minimum_required(VERSION 3.5)

project(gray_scott_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(gray_scott_cpp
    src/main.cpp
    src/linalg.cpp src/linalg.h
    src/solver.cpp src/solver.h
    src/simplecmdparser.hpp
    src/grayscott.h src/grayscott.cpp
)

target_link_libraries(gray_scott_cpp tbb pthread)

install(TARGETS gray_scott_cpp
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
