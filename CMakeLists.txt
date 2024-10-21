cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C CXX)

cmake_policy(SET CMP0167 NEW)

set(CMAKE_CXX_STANDARD 11)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)

# --- Detect SIMD --------------------------------------------------------------

include("src/scripts/cmake/CythonExtension.cmake")
include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindSSE2.cmake")
include("src/scripts/cmake/FindSSSE3.cmake")

# --- Find Boost::math ---------------------------------------------------------

set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(BOOST_MATH_STANDALONE        ON)

find_package(Boost 1.45 COMPONENTS math_c99)

# --- Setup include folders ----------------------------------------------------

include_directories(${CMAKE_CURRENT_LIST_DIR}/vendor/tantan/src)
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

# --- Compile source code ------------------------------------------------------

add_subdirectory(src)
