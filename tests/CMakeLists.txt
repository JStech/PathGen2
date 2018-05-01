find_package(PathGen 0.1 REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
    ${PathGen_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    )

add_executable(wraparound_test wraparound_test.cpp)

target_link_libraries(wraparound_test
    ${PathGen_LIBRARIES}
 )
