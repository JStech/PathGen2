#include <pathgen/IMUGenerator.h>

using namespace pathgen;

void IMUGenerator::GenerateInertialMeasurements(
        const PoseSpline& pose_spline,
        const ImuParameters& imu_parameters,
        const ImuPathOptions& path_options,
        ImuMeasurementDeque* imu_measurements,
        std::vector<PosePtr>* imu_poses)
{
    // generate the imu measurements that correspond to the given trajectory
    int NUM_MEASUREMENTS = (int)(path_options.duration*imu_parameters.rate);
    double spline_increment = 1.0/NUM_MEASUREMENTS; // how many pieces to divide the spline into (spline goes from t: 0->1)
    double spline_pos = 0.0; // initial spline position, t=0;
    const double DT = 1.0/double(imu_parameters.rate);

    double t0 = 0.0;

    for (size_t t = 0; t <= (size_t)NUM_MEASUREMENTS; ++t) {

        // get the current ground truth pose (for debug purposes only)
        imu_poses->push_back(std::make_shared<Pose>(pose_spline(spline_pos)));

        // linear acceleration (no need for velocity as its not measured by the imu)
        PointType accel = pose_spline.position_spline.derivatives
                (spline_pos, 2).col(2);

        // get the current orientation
        Eigen::Quaterniond q_wb = pose_spline.quat_spline(spline_pos);

        Eigen::Vector3d acc = Eigen::Vector3d(accel); // world frame acceleration
        acc = acc / (path_options.duration*path_options.duration); // scale the acceleration so it's compatible with the duration

        acc += Eigen::Vector3d(0,0,imu_parameters.g); // add in gravity, measured in the world frame

        // now rotate acceleration into the body frame (world to body frame)
        acc = q_wb.toRotationMatrix().inverse() * acc;

        if(path_options.add_noise_to_imu_measurements){
            // now add noise to imu accel measurements
            acc += Eigen::Vector3d::Random() * imu_parameters.sigma_a_c /
                    sqrt(DT);
        }

        // get the current angular velocity
        Eigen::Vector3d gyro_w = pose_spline.quat_spline.derivative(spline_pos)
                /path_options.duration;

        // and convert from world to body rates
        Eigen::Vector3d pqr = R2pqr(q_wb.toRotationMatrix());
        Eigen::Matrix3d W(Eigen::Matrix3d::Identity());
        W(0, 2) = -sin(pqr[1]);
        W(1, 1) = cos(pqr[0]); W(1, 2) = cos(pqr[1])*sin(pqr[0]);
        W(2, 1) = -sin(pqr[0]); W(2, 2) = cos(pqr[1])*cos(pqr[0]);
        Eigen::Vector3d gyro_body = W*gyro_w;

        if(path_options.add_noise_to_imu_measurements){
            // add noise to gyro measurements
            gyro_body += Eigen::Vector3d::Random() * imu_parameters.sigma_g_c
                    / sqrt(DT);
        }

        imu_measurements->push_back(
                    ImuMeasurement(t0 + t*DT, ImuSensorReadings(gyro_body, acc)));

        spline_pos += spline_increment;

    }

}

void IMUGenerator::SaveMeasurementFiles(const ImuMeasurementDeque& meas){
    // write out imu measurements to a file
    std::ofstream imu_file("imu_meas.txt", std::ios_base::trunc);
    std::ofstream imu_accel("accel.txt", std::ios_base::trunc);
    std::ofstream imu_gyro("gyro.txt", std::ios_base::trunc);
    std::ofstream imu_ts("timestamp.txt", std::ios_base::trunc);
    for(const auto& m : meas){
        imu_file << static_cast<uint64_t>(1e9*m.timestamp) << ", " <<
          m.measurement.accelerometers[0] << "," <<
          m.measurement.accelerometers[1] << ", " <<
          m.measurement.accelerometers[2] << ", " <<
          m.measurement.gyroscopes[0] << ", " <<
          m.measurement.gyroscopes[1] << ", " <<
          m.measurement.gyroscopes[2] << std::endl;

        // measurement files for reading in with HAL
        imu_accel << m.measurement.accelerometers[0] << "," <<
        m.measurement.accelerometers[1] << ", " <<
        m.measurement.accelerometers[2] << std::endl;

        imu_gyro <<  m.measurement.gyroscopes[0] << ", "
         << m.measurement.gyroscopes[1] << ", " <<
        m.measurement.gyroscopes[2] << std::endl;

        imu_ts <<  m.timestamp  << ", " << m.timestamp << std::endl;
    }

}


