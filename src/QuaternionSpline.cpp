#include <pathgen/QuaternionSpline.h>
#include <pathgen/Utils.h>
#include <pathgen/QuaternionAddons.h>

namespace pathgen {

  // quadilatieral quaternion interpolation
  void squad(const Eigen::Quaterniond& q1,
      const Eigen::Quaterniond& q2,
      const Eigen::Quaterniond& s1,
      const Eigen::Quaterniond& s2,
      double t,
      Eigen::Quaterniond* dst) {
    Eigen::Quaterniond dstQ;
    Eigen::Quaterniond dstS;

    dstQ = q1.slerp(t, q2); // use eigen's slerp method
    dstS = s1.slerp(t, s2);
    *dst = dstQ.slerp(2.0f*t*(1.0f-t), dstS);
  }

  // derivatrive of squat w.r.t. t
  void squad_prime(const Eigen::Quaterniond& q1,
      const Eigen::Quaterniond& q2,
      const Eigen::Quaterniond& s1,
      const Eigen::Quaterniond& s2,
      double t,
      Eigen::Quaterniond* dst) {

    if ((t < -1e-4) || (t > 1 + 1e-4)) {
      std::cerr << "Squad_prime(p,a,b,q, t): " << t << " t < 0 or t > 1." <<
        std::endl;
    }

    Eigen::Quaterniond U = q1.slerp(t, q2),
      V = s1.slerp(t, s2),
      W = U.inverse()*V,
      U_prime = U*log(q1.inverse()*q2),
      V_prime = V*log(s1.inverse()*s2),
      W_prime = minus(U.inverse()*V_prime, power(U, -1)*U_prime*W);

    *dst = plus(U *
        plus(power(W, (2*t*(1-t))) * times(log(W), (2-4*t)), power(W,
            (2*t*(1-t)-1)) * times(W_prime, 2*t*(1-t))), U_prime *
        (power(W, (2*t*(1-t)))));
  }

  // derivatrive of slerp w.r.t. t
  void slerp_prime(const Eigen::Quaterniond& q0,
      const Eigen::Quaterniond& q1,
      double t,
      Eigen::Quaterniond* dst) {
    if ((t < 0) || (t > 1))
      std::cerr << "Slerp_prime(q0, q1, t): t < 0 or t > 1. t is set to 0." <<
        std::endl;

    *dst = q0.slerp(t, q1)*log((q0.inverse()*q1));
  }

  // Return angular velocity from a quaternion and it's time derivative
  Eigen::Vector3d omega(const Eigen::Quaterniond& q,
      const Eigen::Quaterniond& q_dot) {
    return 2*(q_dot*q.conjugate()).vec();
  }

  bool TESTsquad_prime(const Eigen::Quaterniond& q1,
      const Eigen::Quaterniond& q2,
      const Eigen::Quaterniond& s1,
      const Eigen::Quaterniond& s2,
      double t, double eps = 1e-8) {
    Eigen::Quaterniond q_prime, q_prime_finite_diff;
    Eigen::Quaterniond q_plus, q_minus;
    squad(q1, q2, s1, s2, t+eps, &q_plus);
    squad(q1, q2, s1, s2, t-eps, &q_minus);
    q_prime_finite_diff.coeffs() = (q_plus.coeffs() - q_minus.coeffs())/(2*eps);
    squad_prime(q1, q2, s1, s2, t, &q_prime);

    double diff = (q_prime_finite_diff.coeffs() - q_prime.coeffs()).norm();
    if (diff > 1e-2) {
      std::cerr << "squad derivative not close enough: " << diff <<
        " t: " << t << std::endl;
      return false;
    }
    return true;
  }

  QuaternionSpline::QuaternionSpline(
      Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> ctrl,
      Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> support,
      Eigen::VectorXd knots) :
    control_pts_(ctrl),
    support_pts_(support),
    knots_(knots) {}

  void QuaternionSpline::SetIntegrateDerivativeForOrientation(bool integrate) {

    if (integrate && integrated_quaternions_.size() == 0) {
      // need to populate integration map

      const double step_size = eps_;
      const double num_integration_steps = 1.0/step_size;
      integrated_quaternions_.reserve(num_integration_steps+1);

      integrated_quaternions_.push_back(control_pts_(0));

      Eigen::Vector3d euler = R2pqr(control_pts_(0).toRotationMatrix());

      for (double t = step_size; t-1.0 <= 1e-6; t+=step_size) {
        Eigen::Vector3d w = derivative(t); // get the angular rates at point t
        euler += w*step_size; // integrate up to that point
        integrated_quaternions_.push_back(euler2quat(euler));
      }

      std::cout << "integrated quaternion and generated " <<
        integrated_quaternions_.size() << " elements" << std::endl;
    }
    integrate_derivative_for_orientation_ = integrate;
  }

  // the spline value at a given locaiton
  // see: http://web.mit.edu/2.998/www/QuaternionReport1.pdf
  Eigen::Quaterniond QuaternionSpline::operator()(const double& u) const {

    Eigen::Quaterniond interpolated_quaternion;

    if (!integrate_derivative_for_orientation_) {

      size_t index = 0;
      for (size_t ii = 0; ii < (size_t)knots_.size()-1; ++ii) {

        if (u < knots_(ii+1)) {
          index = ii;
          break;
        }

        if (ii == (size_t)knots_.size()-2) {
          index = ii+1;
        }
      }

      double size = knots_(index+1) - knots_(index);
      double t = (u - knots_(index))/size;

      // actually do the interpolation now (squad)
      Eigen::Quaterniond q = control_pts_(index);
      Eigen::Quaterniond q_plus_1 = control_pts_(index+1);
      Eigen::Quaterniond s = support_pts_(index);
      Eigen::Quaterniond s_plus_1 = support_pts_(index+1);

      squad(q, q_plus_1, s, s_plus_1, t, &interpolated_quaternion);
    } else {

      // get the pre-integrated values
      size_t idx = (size_t)round(u/eps_);
      if (idx > integrated_quaternions_.size()) {
        std::cerr << "tried to get quaternion out of bounds, idx: " << idx <<
          std::endl;
        interpolated_quaternion = Eigen::Quaterniond::Identity();
      } else {
        interpolated_quaternion = integrated_quaternions_[idx];
      }
    }

    return interpolated_quaternion;
  }


  // see: https://www.geometrictools.com/Documentation/Quaternions.pdf
  Eigen::Quaterniond QuaternionSpline::quaternionDerivative(
      const double& u, Eigen::Quaterniond* quat_ptr) const {
    size_t index = 0;
    for (size_t ii = 0; ii < (size_t)knots_.size()-1; ++ii) {

      if (u < knots_(ii+1)) {
        index = ii;
        break;
      }

      if (ii == (size_t)knots_.size()-2) {
        index = ii+1;
      }
    }

    double size = knots_(index+1) - knots_(index);
    double t = (u - knots_(index))/size;

    Eigen::Quaterniond q = control_pts_(index);
    Eigen::Quaterniond q_plus_1 = control_pts_(index+1);
    Eigen::Quaterniond s = support_pts_(index);
    Eigen::Quaterniond s_plus_1 = support_pts_(index+1);

    Eigen::Quaterniond quat_dot;
    Eigen::Quaterniond quat;
    if (quat_ptr == nullptr) {
      quat_ptr = &quat;
    }

    // get the interpolated quaternion
    squad(q, q_plus_1, s, s_plus_1, t, quat_ptr);
    // and it's derivative w.r.t. t
    squad_prime(q, q_plus_1, s, s_plus_1, t, &quat_dot);

    TESTsquad_prime(q, q_plus_1, s, s_plus_1, t);

    Eigen::Vector4d coeffs = quat_dot.coeffs();
    coeffs*= (control_pts_.size()-1);
    quat_dot.coeffs() = coeffs;

    return quat_dot;
  }


  // see: https://www.geometrictools.com/Documentation/Quaternions.pdf
  Eigen::Vector3d QuaternionSpline::derivative(const double& u) const {
    Eigen::Quaterniond quat;
    Eigen::Quaterniond quat_dot = quaternionDerivative(u, &quat);
    Eigen::Vector3d w = omega(quat, quat_dot);

    return w;
  }


  // inspired by eigen's spline interpolation class (static)
  QuaternionSpline QuaternionSplineFitting::Interpolate(
      const Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> pts,
      const Eigen::VectorXd& knots) {

    Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> support_pts(pts.cols());

    // compute the support points for each control point
    for (size_t ii = 0; ii < (size_t)pts.cols(); ++ii) {
      Eigen::Quaterniond q_i = pts(ii);

      Eigen::Quaterniond q_i_plus_1 = Eigen::Quaterniond::Identity();
      if (ii == (size_t)pts.cols() - 1) {
        q_i_plus_1 = pts(ii);
      } else {
        q_i_plus_1 = pts(ii+1);
      }

      Eigen::Quaterniond q_i_minus_1 = Eigen::Quaterniond::Identity();
      if (ii == 0) {
        q_i_minus_1 = pts(ii);
      } else {
        q_i_minus_1 = pts(ii-1);
      }

      Eigen::Vector3d log_vec = log(q_i.inverse() * q_i_plus_1).vec() +
        log(q_i.inverse() * q_i_minus_1).vec();
      Eigen::Quaterniond s_i = q_i * exp(-1 * log_vec/4.0);
      support_pts(ii) = s_i;
    }

    return QuaternionSpline(pts, support_pts, knots);
  }

} // namespace pathgen
