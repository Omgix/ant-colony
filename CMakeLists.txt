cmake_minimum_required(VERSION 3.5)
project(ant_colony)

set(CMAKE_CXX_STANDARD 14)

include(${PROJECT_SOURCE_DIR}/cmake/misc.cmake)
include_directories(${PROJECT_SOURCE_DIR}/include)

configure_file(${PROJECT_SOURCE_DIR}/cmake/get_eigen.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/eigen-download/CMakeLists.txt)
setup_eigen()

include(CTest)
add_subdirectory( test )

add_library(AntColonyBase AntColonyBase.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(AntColonyBase PUBLIC OpenMP::OpenMP_CXX Eigen3::Eigen)
else()
    target_link_libraries(AntColonyBase Eigen3::Eigen)
endif()