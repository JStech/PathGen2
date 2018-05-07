#pragma once
#ifndef INCLUDE_PATHGEN_QUATERION_SPLINE_H_
#define INCLUDE_PATHGEN_QUATERION_SPLINE_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace pathgen
{

void squad(const Eigen::Quaterniond& q1,
           const Eigen::Quaterniond& q2,
           const Eigen::Quaterniond& s1,
           const Eigen::Quaterniond& s2,
           double t,
           Eigen::Quaterniond* dst);

void squad_prime(const Eigen::Quaterniond& q1,
                 const Eigen::Quaterniond& q2,
                 const Eigen::Quaterniond& s1,
                 const Eigen::Quaterniond& s2,
                 double t,
                 Eigen::Quaterniond* dst);


class QuaternionSpline
{
public:
    QuaternionSpline(){}
    QuaternionSpline(Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> ctrl,
                     Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> support);

    // the spline value at a given locaiton
    Eigen::Quaterniond operator()(const double& u) const;

    Eigen::Vector3d derivative(const double& u) const;
    Eigen::Quaterniond quaternionDerivative(const double& u,
        Eigen::Quaterniond* = nullptr) const;
    void SetIntegrateDerivativeForOrientation(bool integrate);

protected:
    Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> control_pts_; //q_i
    Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> support_pts_; //s_i
    Eigen::VectorXd knots_;
    std::vector<Eigen::Quaterniond> integrated_quaternions_;
    bool integrate_derivative_for_orientation_;
    double eps_ = 1e-3;
};

struct QuaternionSplineFitting
{
    static void ChordLengths(
            const Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic>& pts,
            Eigen::VectorXd& chord_lengths);

    static QuaternionSpline Interpolate(
            const Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> pts);

};


}

#endif
