#include <pathgen/QuaternionSpline.h>
#include <pathgen/Utils.h>
#include <pathgen/QuaternionAddons.h>
#include <glog/logging.h>
#include <iomanip>

namespace pathgen
{

// quadilatieral quaternion interpolation
void squad(const Eigen::Quaterniond& q1,
           const Eigen::Quaterniond& q2,
           const Eigen::Quaterniond& s1,
           const Eigen::Quaterniond& s2,
           double t,
           Eigen::Quaterniond* dst)
{

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
                 Eigen::Quaterniond* dst)
{

    if( (t < -1e-4) || (t > 1 + 1e-4) ){
        std::cerr << "Squad_prime(p,a,b,q, t): " << t << " t < 0 or t > 1." << std::endl;
    }

    Eigen::Quaterniond U = q1.slerp(t, q2),
            V = s1.slerp(t, s2),
            W = U.inverse()*V,
            U_prime = U*log(q1.inverse()*q2),
            V_prime = V*log(s1.inverse()*s2),
            W_prime = minus(U.inverse()*V_prime, power(U, -1)*U_prime*W);

    *dst = plus(U *
                plus( power(W, (2*t*(1-t))) * times(log(W), (2-4*t)), power(W, (2*t*(1-t)-1)) * times(W_prime, 2*t*(1-t) ) )
                , U_prime * ( power(W,(2*t*(1-t))) ));

}

// derivatrive of slerp w.r.t. t
void slerp_prime(const Eigen::Quaterniond& q0,
                 const Eigen::Quaterniond& q1,
                 double t,
                 Eigen::Quaterniond* dst)
{
    if( (t < 0) || (t > 1) )
        std::cerr << "Slerp_prime(q0, q1, t): t < 0 or t > 1. t is set to 0." << std::endl;

    if(q0.dot(q1) >= 0)
        *dst = q0.slerp(t, q1)*log((q0.inverse()*q1));
    else
        *dst = q0.slerp(t, q1)*log(times(q0.inverse(),-1)*q1);
}

// Return angular velocity from a quaternion and it's time derivative
Eigen::Vector3d omega(const Eigen::Quaterniond& q, const Eigen::Quaterniond& q_dot)
{
        Eigen::Matrix3d A;
        Eigen::Vector3d w, B;
        A = 0.5*E(q, 0);
        B = q_dot.vec();
        if(A.determinant())
        {
            w = A.colPivHouseholderQr().solve(B);
            bool solution_exists = (A*w).isApprox(B, 1e-5);
            if(!solution_exists) {
                std::cerr << "QR solution failed." << std::endl;
                std::cerr <<
                  q.w() << " " << q.x() << " " << q.y() << " " << q.z() << ", " <<
                  q_dot.w() << " " << q_dot.x() << " " << q_dot.y() << " " << q_dot.z() << ", " <<
                  std::endl << A << std::endl << B << std::endl << w;
            }
        }
        else
            w.setZero();

//    Eigen::Vector3d w(Eigen::Vector3d::Zero());
//    Eigen::Quaterniond q_w = q_dot*q.inverse();
//    w = q_w.vec()*2;

    return w;
}

bool TESTsquad_prime(const Eigen::Quaterniond& q1,
                     const Eigen::Quaterniond& q2,
                     const Eigen::Quaterniond& s1,
                     const Eigen::Quaterniond& s2,
                     double t, double eps = 1e-8)
{
    Eigen::Quaterniond q_prime, q_prime_finite_diff;
    Eigen::Quaterniond q_plus, q_minus;
    squad(q1, q2, s1, s2, t+eps, &q_plus);
    squad(q1, q2, s1, s2, t-eps, &q_minus);
    q_prime_finite_diff.coeffs() = (q_plus.coeffs() - q_minus.coeffs())/(2*eps);
    squad_prime(q1, q2, s1, s2, t, &q_prime);

    double diff = (q_prime_finite_diff.coeffs() - q_prime.coeffs()).norm();
    if(diff > 1e-2){
        std::cerr << "squad derivative not close enough: " << diff <<
                     " t: " << t << std::endl;
        return false;
    }
    return true;
}

QuaternionSpline::QuaternionSpline(Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> ctrl,
                                   Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> support):
    control_pts_(ctrl),
    support_pts_(support)
{
    // build the chord lengths
   QuaternionSplineFitting::ChordLengths(ctrl, knots_);
}


void QuaternionSpline::SetIntegrateDerivativeForOrientation(bool integrate){

    if(integrate && integrated_quaternions_.size() == 0){
        // need to populate integration map

        const double step_size = eps_;
        const double num_integration_steps = 1.0/step_size;
        integrated_quaternions_.reserve(num_integration_steps+1);

        integrated_quaternions_.push_back(control_pts_(0));

        Eigen::Vector3d euler = R2pqr(control_pts_(0).toRotationMatrix());

        for(double t = step_size; t-1.0<=1e-6; t+=step_size){
            Eigen::Vector3d w = derivative(t); // get the angular rates at point t
            euler += w*step_size; // integrate up to that point
            integrated_quaternions_.push_back(euler2quat(euler));
        }

        std::cout << "integrated quaternion and generated " << integrated_quaternions_.size() << " elements" << std::endl;

    }
    integrate_derivative_for_orientation_ = integrate;
}

// the spline value at a given locaiton
//see: http://web.mit.edu/2.998/www/QuaternionReport1.pdf
Eigen::Quaterniond QuaternionSpline::operator()(const double& u) const
{

    Eigen::Quaterniond interpolated_quaternion;

    if(!integrate_derivative_for_orientation_){

        size_t index = 0;
        for(size_t ii = 0; ii < (size_t)knots_.size()-1; ++ii){

            if(u < knots_(ii+1)){
                index = ii;
                break;
            }

            if(ii == (size_t)knots_.size()-2)
                index = ii+1;
        }


        double size = knots_(index+1) - knots_(index);
        double t = (u - knots_(index))/size;

        // actually do the interpolation now (squad)
        Eigen::Quaterniond q = control_pts_(index);
        Eigen::Quaterniond q_plus_1 = control_pts_(index+1);
        Eigen::Quaterniond s = support_pts_(index);
        Eigen::Quaterniond s_plus_1 = support_pts_(index+1);

        squad(q, q_plus_1, s, s_plus_1, t, &interpolated_quaternion);

    }else{
        // this is hacky...integrate the 1st order derivative to get the orientation.
        // needed due to inconsistent derivative and spline for trajectories with small angular velocity
        // TODO: fix the derivatives so this isn't needed.

        /*
        Eigen::Vector3d euler = R2pqr(control_pts_(0).toRotationMatrix());
        const double step_size = eps_;

        for(double t = step_size; t<=u; t+=step_size){
            Eigen::Vector3d w = derivative(t);
            euler += w*step_size;
        }

        interpolated_quaternion = euler2quat(euler);
        */

        // If integrating quaternion_dot directly, then it behaves as expected and yields
        // the original curve, something is wrong in getting the angular rates out of the quaternion
        // derivative...TODO (q_dot = 1/2 * w * q)

        /*
        Eigen::Quaterniond q = control_pts_(0);
        const double step_size = 1e-3;
        for(double t = step_size; t <=u; t+=step_size){
            q.coeffs() = q.coeffs() + step_size*quaternionDerivative(t).coeffs();
            q.normalize();
        }

        interpolated_quaternion = q;
        */

        // get the pre-integrated values
        size_t idx = (size_t)round(u/eps_);
        if( idx > integrated_quaternions_.size()){
            std::cerr << "tried to get quaternion out of bounds, idx: " << idx << std::endl;
            interpolated_quaternion = Eigen::Quaterniond::Identity();
        }else{
            interpolated_quaternion = integrated_quaternions_[idx];
        }

    }

    return interpolated_quaternion;

}


// see: https://www.geometrictools.com/Documentation/Quaternions.pdf
Eigen::Quaterniond QuaternionSpline::quaternionDerivative(const double& u) const
{
    size_t index = 0;
    for(size_t ii = 0; ii < (size_t)knots_.size()-1; ++ii){

        if(u < knots_(ii+1)){
            index = ii;
            break;
        }

        if(ii == (size_t)knots_.size()-2)
            index = ii+1;
    }

    double size = knots_(index+1) - knots_(index);
    double t = (u - knots_(index))/size;

    Eigen::Quaterniond q = control_pts_(index);
    Eigen::Quaterniond q_plus_1 = control_pts_(index+1);
    Eigen::Quaterniond s = support_pts_(index);
    Eigen::Quaterniond s_plus_1 = support_pts_(index+1);

    Eigen::Quaterniond quat_dot;
    Eigen::Quaterniond quat;


    squad(q, q_plus_1, s, s_plus_1, t, &quat); // get the interpolated quaternion
    squad_prime(q, q_plus_1, s, s_plus_1, t, &quat_dot); // and it's derivative w.r.t. t

    Eigen::Vector4d coeffs = quat_dot.coeffs();
    coeffs*= (control_pts_.size()-1);
    quat_dot.coeffs() = coeffs;

    return quat_dot;
}


// see: https://www.geometrictools.com/Documentation/Quaternions.pdf
Eigen::Vector3d QuaternionSpline::derivative(const double& u) const
{
    size_t index = 0;
    for(size_t ii = 0; ii < (size_t)knots_.size()-1; ++ii){

        if(u < knots_(ii+1)){
            index = ii;
            break;
        }

        if(ii == (size_t)knots_.size()-2)
            index = ii+1;
    }

    double size = knots_(index+1) - knots_(index);
    double t = (u - knots_(index))/size;

    Eigen::Quaterniond q = control_pts_(index);
    Eigen::Quaterniond q_plus_1 = control_pts_(index+1);
    Eigen::Quaterniond s = support_pts_(index);
    Eigen::Quaterniond s_plus_1 = support_pts_(index+1);

    Eigen::Quaterniond quat_dot;
    Eigen::Quaterniond quat;

    //    if(index == 1 || index == knots_.size()-1 ){
    //        // we're between the first and second keyframe...do slerp since squad
    //        // needs more points
    //        // or we're between  next-to-last and last keyframe
    //        quat = q.slerp(t, q_plus_1);
    //        slerp_prime(q, q_plus_1, t, &quat_dot);
    //    }else {
    //        // we're somewhere in the middle, can do squad
    //        squad(q, q_plus_1, s, s_plus_1, t, &quat);
    //        squad_prime(q, q_plus_1, s, s_plus_1, t, &quat_dot);
    //    }

    squad(q, q_plus_1, s, s_plus_1, t, &quat); // get the interpolated quaternion
    squad_prime(q, q_plus_1, s, s_plus_1, t, &quat_dot); // and it's derivative w.r.t. t

    // LOG(INFO) << std::fixed << std::setprecision(8) << u << " " <<
    //   q.w() << " " << q.x() << " " << q.y() << " " << q.z() << ", " <<
    //   q_plus_1.w() << " " << q_plus_1.x() << " " << q_plus_1.y() << " " << q_plus_1.z() << ", " <<
    //   s.w() << " " << s.x() << " " << s.y() << " " << s.z() << ", " <<
    //   s_plus_1.w() << " " << s_plus_1.x() << " " << s_plus_1.y() << " " << s_plus_1.z() << ", " <<
    //   quat_dot.w() << " " << quat_dot.x() << " " << quat_dot.y() << " " << quat_dot.z() << ", " <<
    //   quat.w() << " " << quat.x() << " " << quat.y() << " " << quat.z();

    TESTsquad_prime(q, q_plus_1, s, s_plus_1, t);

    //quat_dot.vec() *= (control_pts_.size()-1);
    Eigen::Vector4d coeffs = quat_dot.coeffs();
    coeffs*= (control_pts_.size()-1);
    quat_dot.coeffs() = coeffs;
    Eigen::Vector3d w = omega(quat, quat_dot);

    //sanity check for w: q_dot = 1/2 * w * q;
    //    double diff = (quat_dot.coeffs() - 0.5*(Eigen::Quaterniond(0, w[0], w[1], w[2])*quat).coeffs()).norm();
    //    if(diff > 1e-3){
    //        std::cerr << "omega not consistent, diff: " << diff << std::endl;
    //    }

    return w;


    //    // do finite differences for now
    //    Eigen::Quaterniond dq;
    //    Eigen::Vector3d w;
    //    const double eps = 1e-8;
    //    if(t >= eps && t <= 1-eps){
    //        Eigen::Quaterniond qp;
    //        Eigen::Quaterniond qm;
    //        squad(q, q_plus_1, s, s_plus_1, (t+eps), &qp);
    //        squad(q, q_plus_1, s, s_plus_1, (t-eps), &qm);
    //        Eigen::Vector3d euler_plus = qp.toRotationMatrix().eulerAngles(0,1,2);
    //        Eigen::Vector3d euler_minus = qm.toRotationMatrix().eulerAngles(0,1,2);
    //        w = (euler_plus- euler_minus)/(2*eps);
    //    }else if (t < eps){
    //        Eigen::Vector3d temp_vec = log(q.inverse()*q_plus_1).vec() + 2*log(q.inverse()*s).vec();
    //        Eigen::Quaterniond temp(0, temp_vec[0], temp_vec[1], temp_vec[2] );
    //        dq = q*temp;
    //        w = dq.toRotationMatrix().eulerAngles(0,1,2);
    //    }else{
    //        Eigen::Vector3d temp_vec = log(q.inverse()*q_plus_1).vec() - 2*log(q_plus_1.inverse()*s_plus_1).vec();
    //        Eigen::Quaterniond temp(0, temp_vec[0], temp_vec[1], temp_vec[2] );
    //        dq = q_plus_1*temp;
    //        w = dq.toRotationMatrix().eulerAngles(0,1,2);
    //    }

    //    //    Eigen::Quaterniond U = q.slerp(t, q_plus_1);
    //    //    Eigen::Quaterniond V = s.slerp(t, s_plus_1);
    //    //    Eigen::Quaterniond dU= U*log(q.inverse()*q_plus_1);
    //    //    Eigen::Quaterniond dV = V*log(s.inverse()*s_plus_1);
    //    //    Eigen::Quaterniond W = U.inverse()*V;
    //    //    Eigen::Quaterniond Wexp = dU * exp( log(W) * (2*t*(1-t)) );

    //    //    double a = acos(W.w());
    //    //    double sina = sin(a);
    //    //    Eigen::Vector3d v;
    //    //    v.setZero();
    //    //    if (sina > 0)
    //    //    {
    //    //        v[0] = W.x()/sina;
    //    //        v[1] = W.y()/sina;
    //    //        v[2] = W.z()/sina;
    //    //    }



    //    return w;
}


// inspired by eigen's spline interpolation class (static)
QuaternionSpline QuaternionSplineFitting::Interpolate(const Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> pts)
{

    Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> support_pts(pts.cols());

    // compute the support points for each control point
    for(size_t ii = 0; ii < (size_t)pts.cols(); ++ii)
    {
        Eigen::Quaterniond q_i = pts(ii);

        Eigen::Quaterniond q_i_plus_1 = Eigen::Quaterniond::Identity();
        if(ii == (size_t)pts.cols() - 1)
            q_i_plus_1 = pts(ii);
        else
            q_i_plus_1 = pts(ii+1);

        Eigen::Quaterniond q_i_minus_1 = Eigen::Quaterniond::Identity();
        if(ii == 0)
            q_i_minus_1 = pts(ii);
        else
            q_i_minus_1 = pts(ii-1);

        Eigen::Vector3d log_vec = log(q_i.inverse() * q_i_plus_1).vec() + log(q_i.inverse() * q_i_minus_1).vec();
        Eigen::Quaterniond s_i = q_i * exp(-1 * log_vec/4.0);
        support_pts(ii) = s_i;
        // LOG(INFO) << std::fixed << std::setprecision(8) <<
        //   " " << pts(ii).w() << " " << pts(ii).x() <<
        //   " " << pts(ii).y() << " " << pts(ii).z() << "," <<
        //   " " << support_pts(ii).w() << " " << support_pts(ii).x() <<
        //   " " << support_pts(ii).y() << " " << support_pts(ii).z();

    }

    return QuaternionSpline(pts, support_pts);


}

// inspired by eigen's spline interpolation class
void QuaternionSplineFitting::ChordLengths(const Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic>& pts,
                                          Eigen::VectorXd& chord_lengths)
{
    const int n = pts.cols();

    chord_lengths.resize(pts.cols());
    chord_lengths[0] = 0;

    for(size_t i = 0; i < (size_t)n-1 ; ++i){
        // get the angle between quaternions
        chord_lengths(i+1) = pts(i).dot(pts(i+1));
        if (chord_lengths(i+1) < 0) {
          LOG(INFO) << "Fixing negative chord length, " << i+1 << " " <<
            chord_lengths(i+1);
          chord_lengths(i+1) *= -1;
        }
    }

    std::partial_sum(chord_lengths.data(), chord_lengths.data()+n, chord_lengths.data());

    chord_lengths /= chord_lengths(n-1); // normalize
}

}
