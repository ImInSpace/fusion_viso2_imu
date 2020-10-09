cmake_minimum_required(VERSION 3.14)
project(fusion_viso2_imu)

set(USERNAME $ENV{USERNAME})
message(${USERNAME})
if(${USERNAME} STREQUAL "Inigo Moreno")
set(USERNAME "Moreno Inigo")
endif()
message(${USERNAME})

#include eigen
include_directories("C:\\Users\\${USERNAME}\\Documents\\libraries\\eigen-3.3.7")

#include boost
set(BOOST_ROOT "C:\\Users\\${USERNAME}\\Documents\\libraries\\boost_1_60_0")
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

#include yaml-cpp
add_subdirectory(third-party/yaml-cpp)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

#include library
add_subdirectory(library)
add_subdirectory(functions)
add_executable(test_fusion_viso2_imu test.cpp multivar_noise.h)
add_executable(test_rotations test_rotations.cpp)

#link libraries
target_link_libraries(test_fusion_viso2_imu PUBLIC fusion_viso2_imu functions_fusion_viso2_imu)
target_include_directories(test_fusion_viso2_imu PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(test_fusion_viso2_imu ${Boost_LIBRARIES})
target_link_libraries(test_fusion_viso2_imu PRIVATE yaml-cpp)