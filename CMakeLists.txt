cmake_minimum_required(VERSION 3.16)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Use local version of TinyAD
add_subdirectory(${PROJECT_SOURCE_DIR}/TinyAD)
message(STATUS "tinyad_SOURCE_DIR: ${PROJECT_SOURCE_DIR}/TinyAD")

# include(FetchContent)
# FetchContent_Declare(
#   tinyad
#   GIT_REPOSITORY https://github.com/patr-schm/tinyad.git
#   GIT_TAG 75093e14ef0d7bb39657c5f3b2aba1251afaa38c
# )
# #FetchContent_MakeAvailable(tinyad)

# FetchContent_GetProperties(tinyad)
# if(NOT tinyad_POPULATED)
#   # Fetch the content using previously declared details
#   FetchContent_Populate(tinyad)
#   message(STATUS "tinyad_SOURCE_DIR: ${tinyad_SOURCE_DIR}")
#   message(STATUS "tinyad_BINARY_DIR: ${tinyad_BINARY_DIR}")
#   add_subdirectory(${tinyad_SOURCE_DIR} ${tinyad_BINARY_DIR})
# endif()

#include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

#include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Add your project files
file(GLOB SRC_FILES main.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES} ${H_FILES} ${CPP_FILES} ${H_FILES_LBFGS})
target_compile_definitions(${PROJECT_NAME} PUBLIC SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw TinyAD)
# tinyad needs C++17 (not sure if that means this project does too)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
