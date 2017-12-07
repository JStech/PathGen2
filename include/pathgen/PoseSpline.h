#pragma once
#ifndef INCLUDE_PATHGEN_POSE_SPLINE_H_
#define INCLUDE_PATHGEN_POSE_SPLINE_H_

#include <Eigen/Core>
#include "QuaternionSpline.h"
#include "Types.h"
#include "sophus/se3.hpp"
#include <iostream>

namespace pathgen{
struct PoseSpline
{
    PoseSpline(){}

    void SetQuatSpline(QuaternionSpline quat_spline_new){
        quat_spline = quat_spline_new;
    }

    void SetPositionSpline(Eigen::Spline3d position_spline_new){
        position_spline = position_spline_new;
    }

    PoseSpline(QuaternionSpline quat_spline,
               Eigen::Spline3d position_spline):
        quat_spline(quat_spline),
        position_spline(position_spline){}

    Pose operator()(double u) const
    {

        Eigen::Quaterniond q = quat_spline(u);
        Eigen::Vector3d pos = position_spline(u);

        // now build transformation with rotation and position
        return Pose(q, pos);
    }

    pathgen::QuaternionSpline quat_spline;
    Eigen::Spline3d position_spline;
};
}

#endif
