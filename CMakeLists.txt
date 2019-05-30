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
target_link_libraries(AntColonyBase Eigen3::Eigen)

SET( PACKAGE_VERSION "0.1.1" )
FIND_PACKAGE( pybind11 REQUIRED )
FIND_PACKAGE(PythonInterp 3)
FIND_PACKAGE(PythonLibs   3)


pybind11_add_module(antcolony SHARED binding.cpp AntColonyBase.cpp)
target_link_libraries(antcolony PRIVATE Eigen3::Eigen)

# The library must not have any prefix and should be located in
# a subfolder that includes the package name. The setup will be
# more complicated otherwise.
SET_TARGET_PROPERTIES( antcolony
        PROPERTIES
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/antcolony"
        SUFFIX ".so"
        )

# This is required for linking the library under Mac OS X. Moreover,
# the suffix ensures that the module can be found by the interpreter
# later on.
IF( APPLE )
    SET_TARGET_PROPERTIES( antcolony
            PROPERTIES
            LINK_FLAGS "-undefined dynamic_lookup"
            SUFFIX     ".so"
            )
ENDIF()


# Place the initialization file in the output directory for the Python
# bindings. This will simplify the installation.
CONFIGURE_FILE( antcolony/__init__.py
        ${CMAKE_CURRENT_BINARY_DIR}/antcolony/__init__.py
        )

# Ditto for the setup file.
CONFIGURE_FILE( antcolony/setup.py
        ${CMAKE_CURRENT_BINARY_DIR}/antcolony/setup.py
        )
