project(mixmat)
cmake_minimum_required(VERSION 2.8)

#add the customised package searches to the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)
find_package(WignerSymbols REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${WignerSymbols_INCLUDE_DIRS})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra")
add_executable(mixmat mixmat.cpp)
target_link_libraries(mixmat ${WignerSymbols_LIBRARIES})
