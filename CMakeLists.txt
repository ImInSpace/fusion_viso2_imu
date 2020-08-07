cmake_minimum_required(VERSION 3.17)
project(fusion_viso2_imu)

include_directories("C:\\Users\\inigo\\Documents\\libraries\\eigen-3.3.7")

set(BOOST_ROOT "C:\\Users\\inigo\\Documents\\libraries\\boost_1_73_0")
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)

add_library(fusion_viso2_imu library.cpp library.h)
add_executable(test_fusion_viso2_imu test.cpp library.cpp library.h multivar_noise.h)