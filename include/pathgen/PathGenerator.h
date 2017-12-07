#pragma once
#ifndef PATHGEN_PATH_GENERATOR_H
#define PATHGEN_PATH_GENERATOR_H

#include <pathgen/Types.h>
#include <pathgen/PoseSpline.h>

namespace pathgen{

class PathGenerator
{
    PathGenerator(){}

public:

    static std::vector<PosePtr> RandomWalk(size_t num_poses,
                                           double radius,
                                           double distance,
                                           double max_speed = 0.1);

    static PoseSpline GenerateSplineFromPoses(std::vector<PosePtr> poses,
                                        bool integrate_gyro = true);
};

}
#endif
