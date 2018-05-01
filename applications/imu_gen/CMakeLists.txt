find_package(Gflags REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Compass REQUIRED)
find_package(PathGen 0.1 REQUIRED)

include_directories(
    ${PathGen_INCLUDE_DIRS}
    ${Compass_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    )

set (SRC
    imu_gen.cpp
    )

add_executable(imu_gen ${SRC})

target_link_libraries(imu_gen
    ${Compass_LIBRARIES}
    ${PathGen_LIBRARIES}
 )
