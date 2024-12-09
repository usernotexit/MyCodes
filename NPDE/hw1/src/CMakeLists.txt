cmake_minimum_required(VERSION 3.25)

project(NPDE1 CXX)

# add_library(Eigen3)
# find_package(Eigen3 REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/../lib/eigen-3.4.0)

add_executable(MyExe main.cpp)
