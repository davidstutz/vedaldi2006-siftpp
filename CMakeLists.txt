cmake_minimum_required(VERSION 2.8)

add_library(siftpp sift.cpp)
add_executable(sift sift-driver.cpp)
target_link_libraries(sift siftpp)
