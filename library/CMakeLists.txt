add_library(fusion_viso2_imu library.cpp library.h)
target_include_directories(fusion_viso2_imu PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(fusion_viso2_imu ${Boost_LIBRARIES})