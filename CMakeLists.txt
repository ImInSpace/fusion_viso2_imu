cmake_minimum_required(VERSION 3.14)
project(fusion_viso2_imu)

include_directories("C:\\Users\\inigo\\Documents\\libraries\\eigen-3.3.7")

set(BOOST_ROOT "C:\\Users\\inigo\\Documents\\libraries\\boost_1_74_0")
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

add_subdirectory(library)
add_executable(test_fusion_viso2_imu test.cpp multivar_noise.h)

target_link_libraries(test_fusion_viso2_imu PUBLIC fusion_viso2_imu)
target_include_directories(test_fusion_viso2_imu PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(test_fusion_viso2_imu ${Boost_LIBRARIES})