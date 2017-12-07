find_package(PathGen 0.1 REQUIRED)
find_package(Pangolin REQUIRED)
find_package(SceneGraph REQUIRED)

include_directories(
    ${PathGen_INCLUDE_DIRS}
    ${Pangolin_INCLUDE_DIRS}
    ${SceneGraph_INCLUDE_DIRS}
    )

set (SRC
    visualizer.cpp
    )

add_executable(visualizer ${SRC})
target_link_libraries(visualizer
    ${PathGen_LIBRARIES}
    ${Pangolin_LIBRARIES}
    ${SceneGraph_LIBRARIES}
 )
