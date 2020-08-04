cmake_minimum_required(VERSION 3.17)
project(fusion_viso2_imu)

include_directories("C:\\Users\\inigo\\Documents\\libraries\\eigen-3.3.7")
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 20)

add_library(fusion_viso2_imu library.cpp library.h)
add_executable(test_fusion_viso2_imu test.cpp library.cpp library.h multivar_noise.h)