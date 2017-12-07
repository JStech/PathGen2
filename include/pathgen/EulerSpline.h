#pragma once
#ifndef INCLUDE_PATHGEN_EULER_SPLINE_H_
#define INCLUDE_PATHGEN_EULER_SPLINE_H_

#include <pathgen/QuaternionSpline.h>
#include <pathgen/Types.h>

namespace pathgen {
struct EulerSpline
{
    EulerSpline(Eigen::ArrayXXd euler_control_pts, int order):
        euler_control_pts_(euler_control_pts){

        for(size_t ii = 0; ii < 3; ++ii){
            Spline1d::ParameterVectorType knots;
            Eigen::ChordLengths(euler_control_pts_.row(ii), knots);
            euler_knots_(ii) = knots;
            Spline1d spline =
                    Eigen::SplineFitting<Spline1d>::Interpolate(euler_control_pts_.row(ii), order);
            splines_(ii) = spline;
        }

    }

    Eigen::Vector3d operator()(const double t){
        Eigen::Vector3d interpolated_pt(Eigen::Vector3d::Zero());

        interpolated_pt[0] = splines_[0](t)(0);
        interpolated_pt[1] = splines_[1](t)(0);
        interpolated_pt[2] = splines_[2](t)(0);


        return interpolated_pt;
    }

    Eigen::MatrixXd derivatives(const double t, int order){
        Eigen::MatrixXd derivatives(3, order+1);
        derivatives.setZero();

        // only filling out the 1st order derivative as this should be all
        // we care about

        derivatives(0, 1) = splines_[0].derivatives(t, order)(0, 1);
        derivatives(1, 1) = splines_[1].derivatives(t, order)(0, 1);
        derivatives(2, 1) = splines_[2].derivatives(t, order)(0, 1);

        return derivatives;
    }

    Eigen::ArrayXXd euler_control_pts_;
    Eigen::Matrix<Spline1d::ParameterVectorType, 3, 1> euler_knots_;
    Eigen::Matrix<Spline1d, 3, 1> splines_;
};
}

#endif
