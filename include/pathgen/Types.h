#pragma once
#ifndef INCLUDE_PATHGEN_TYPES_H_
#define INCLUDE_PATHGEN_TYPES_H_

#include <unsupported/Eigen/Splines>
#include <unsupported/Eigen/AutoDiff>
#include <sophus/se3.hpp>
namespace pathgen
{
typedef Eigen::Spline3d::PointType PointType;
typedef Eigen::Spline3d::KnotVectorType KnotVectorType;
typedef Eigen::Spline<double,1> Spline1d;
typedef Sophus::SE3d Pose;
typedef std::shared_ptr<Sophus::SE3d> PosePtr;

}

#endif /* INCLUDE_PATHGEN_TYPES_H_ */
