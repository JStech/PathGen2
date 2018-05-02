#pragma once
#ifndef INCLUDE_PATHGEN_QUATERION_ADDONS_H_
#define INCLUDE_PATHGEN_QUATERION_ADDONS_H_

#include <Eigen/Core>

namespace pathgen {

  Eigen::Matrix3d x_prod_matrix(const Eigen::Vector3d & x) {
    //!  @brief Cross product matrix for vector x
    Eigen::Matrix3d S;
    S.setZero();
    S(1, 2) = -x(3); S(1, 3) =  x(2);
    S(2, 1) =  x(3); S(2, 3) = -x(1);
    S(3, 1) = -x(2); S(3, 2) =  x(1);

    return S;
  }

  Eigen::Quaterniond log(const Eigen::Quaterniond& q1) {
    Eigen::Quaterniond q;
    q.w() = 0;
    double theta = std::acos(q1.w());
    double sin_theta = std::sin(theta);

    if (std::fabs(sin_theta) > 1e-16) {
      q.vec() = q1.vec()/sin_theta*theta;
    } else {
      q.vec() = q1.vec();
    }

    return q;
  }

  Eigen::Quaterniond exp(const Eigen::Quaterniond& q1) {
    Eigen::Quaterniond q;
    double theta = std::sqrt(q1.vec().dot(q1.vec()));
    double sin_theta = sin(theta);

    q.w() = cos(theta);
    if (std::fabs(sin_theta) > 1e-16) {
      q.vec() = q1.vec()*sin_theta/theta;
    } else {
      q.vec() = q1.vec();
    }

    return q;
  }

  Eigen::Quaterniond exp(const Eigen::Vector3d& v) {
    Eigen::Quaterniond q;
    q.vec() = v;
    q.w() = 0;
    return exp(q);
  }

  Eigen::Quaterniond plus(const Eigen::Quaterniond& q1,
      const Eigen::Quaterniond& q2) {
    Eigen::Quaterniond quat;
    quat.coeffs() = q1.coeffs() + q2.coeffs();
    return quat;
  }

  Eigen::Quaterniond minus(const Eigen::Quaterniond& q1,
      const Eigen::Quaterniond& q2) {
    Eigen::Quaterniond quat;
    quat.coeffs() = q1.coeffs() - q2.coeffs();
    return quat;
  }

  Eigen::Quaterniond times(const Eigen::Quaterniond& q1, const double& c) {
    Eigen::Quaterniond quat;
    quat.coeffs() = q1.coeffs()*c;
    return quat;
  }

  Eigen::Quaterniond power(const Eigen::Quaterniond& q1, const double& t) {
    return exp(times(log(q1), t));
  }

  Eigen::Matrix3d E(const Eigen::Quaterniond& q1, const int16_t sign) {
    Eigen::Matrix3d E, I;
    I.setIdentity();

    if (sign == 1) {
      E = q1.w()*I + x_prod_matrix(q1.vec());
    } else {
      E = q1.w()*I - x_prod_matrix(q1.vec());
    }

    return E;
  }
} // namespace pathgen

#endif
