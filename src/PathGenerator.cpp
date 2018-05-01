#include <pathgen/PathGenerator.h>
#include <glog/logging.h>
#include <pathgen/Utils.h>

using namespace pathgen;

std::vector<PosePtr> PathGenerator::RandomWalk(size_t num_poses,
                                               double radius,
                                               double distance,
                                               double max_speed){

    std::vector<PosePtr> poses;
    poses.push_back(std::make_shared<Pose>()); // first pose is identity;
    Eigen::Vector3d velocity;
    velocity.setZero();


    for(size_t ii = 1; ii < num_poses; ++ii)
    {
        double phi = M_PI * std::fabs(Eigen::Vector2d::Random()[0]);
        double theta = 2*M_PI*std::fabs(Eigen::Vector2d::Random()[0]);

        //        phi = M_PI/2;
        //        theta = -M_PI/9;

        //                phi = (M_PI/2 + M_PI/9);
        //                theta = 0.0;

        double z = radius*cos(phi); // height
        double r = radius*sin(phi); // radius of disk
        double x = r * cos(theta);
        double y = r * sin(theta);

        Eigen::Vector3d forward = poses.back()->rotationMatrix().block<3, 1>(0, 0);
        Eigen::Vector3d position = poses.back()->translation();


        Eigen::Vector3d sphere_pt = Eigen::Vector3d(x, y, z);
        if(poses.size() > 0){
            sphere_pt = poses.back()->rotationMatrix()*sphere_pt;
        }
        Eigen::Vector3d sphere_position = distance * forward + position;
        Eigen::Vector3d sample_position = sphere_position + sphere_pt;
        Eigen::Vector3d force_vec = sample_position - position; // vector from the current position to the point on the sphere we're sampling

        velocity += force_vec; // add to current velocity vector (simulate's inertia)

        if(velocity.norm() > max_speed){
            velocity = (max_speed / velocity.norm()) * velocity;
        }

        // update pose
        forward = velocity.normalized(); // get the direction we want to go in

        Eigen::Vector3d down = -1 * poses.back()->rotationMatrix().block<3, 1>(0, 1).cross(forward);
        down.normalize();
        Eigen::Vector3d right = -1 * forward.cross(down);

        Eigen::Matrix4d new_T;
        new_T.setIdentity();

        new_T.block<3, 1>(0,0) = forward;
        new_T.block<3, 1>(0,1) = right;
        new_T.block<3, 1>(0,2) = down;
        new_T.block<3, 1>(0,3) = poses.back()->translation() + velocity;


        PosePtr new_pose(new Pose(new_T));

        poses.push_back(new_pose);

    }

    return poses;
}

PoseSpline PathGenerator::GenerateSplineFromPoses(std::vector<PosePtr> poses,
                                            bool integrate_gyro){

    // run a cubic spline through the position (for generating imu info)
    Eigen::ArrayXXd positions(3, poses.size());
    Eigen::Array<Eigen::Quaterniond, 1, Eigen::Dynamic> rotations(poses.size());
    Eigen::Quaterniond neg1(-1, 0, 0, 0);
    for(size_t ii = 0; ii < poses.size(); ++ii){
        // translation (x, y, x)
        positions.col(ii) = poses.at(ii)->translation();
        // quaternion (i, j, k, w), eigen ordering
        rotations(ii) = poses.at(ii)->unit_quaternion();
        // ensure our interpolations go the short way around the great arc
        if (ii > 0 && rotations(ii).dot(rotations(ii-1)) < 0) {
          rotations(ii) *= neg1;
        }
    }

    // do the quaternion interpolation
    Eigen::VectorXd quat_chord_lengths;
    quat_chord_lengths.setZero();
    QuaternionSplineFitting::ChordLengths(rotations, quat_chord_lengths);

    QuaternionSpline quat_spline =
            QuaternionSplineFitting::Interpolate(rotations);

    // integrate derivative back so that curves are consistent
    quat_spline.SetIntegrateDerivativeForOrientation(integrate_gyro);


    // And the position interpolation
    // gets the chord lengths corresponding to the points the spline passes through
    Eigen::Spline3d::ParameterVectorType knots;
    Eigen::ChordLengths(positions, knots);

    // 4th order spline so we have smooth accelerations...not really necessary
    Eigen::Spline3d position_spline =
            Eigen::SplineFitting<Eigen::Spline3d>::Interpolate(positions,
                                                               4);

    // Check that the splines are going through the control points
    for(size_t i = 0; i < poses.size(); ++i)
    {
        // check that the quaternion interpolated correctly
        if(!integrate_gyro){
            Eigen::Quaterniond ref_quat = poses.at(i)->unit_quaternion();
            Eigen::Quaterniond spline_quat = quat_spline(quat_chord_lengths(i));
            Eigen::Quaterniond check_quat = ref_quat.inverse()*spline_quat;
            Eigen::Quaterniond neg1(-1, 0, 0, 0);
            if( !check_quat.isApprox(Eigen::Quaterniond::Identity()) &&
                !check_quat.isApprox(neg1)){
                LOG(ERROR) << "Quaternion spline interpolation too far off: " <<
                  check_quat.w() << " " << check_quat.x() << " " <<
                  check_quat.y() << " " << check_quat.z();
            }
        }

        // check that the translation interpolated correctly
        Eigen::Vector3d ref_position = poses.at(i)->translation();
        Eigen::Vector3d spline_position = position_spline(knots(i));
        if( std::fabs((ref_position - spline_position).norm()) > 1e-4 ){
            LOG(ERROR) << "Position spline interpolation too far off: " <<
                          (ref_position - spline_position).transpose();
        }

    }

    // create a pose spline which outputs the full (translatinon + rotation) pose
    return PoseSpline(quat_spline, position_spline);

}

