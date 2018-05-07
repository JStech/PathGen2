#include <pathgen/PathGen.h>
#include <Eigen/Geometry>
#include <glog/logging.h>
#define _USE_MATH_DEFINES
#include <math.h>

uint64_t genPoses(std::vector<pathgen::PosePtr>* poses) {
  Eigen::Matrix4d T;
  Eigen::Quaterniond q(1, 0, 0, 0);

  for (int32_t t = 0; t < 11; t++) {
    double theta = static_cast<double>(t-5) * M_PI / 30;
    T.topRightCorner(3, 1) = Eigen::Vector3d(
        sin(M_PI + theta), cos(M_PI + theta), 0);
    q = Eigen::AngleAxisd(-theta, Eigen::Vector3d::UnitZ());
    T.topLeftCorner(3, 3) = q.toRotationMatrix();
    pathgen::PosePtr new_pose(new pathgen::Pose(T));
    poses->push_back(new_pose);
  }
  return 2000000000;
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  FLAGS_colorlogtostderr = true;

  pathgen::ImuParameters imu_params;

  imu_params.a_max = 176.;
  imu_params.g_max = 7.8;
  imu_params.sigma_g_c = 12e-4;
  imu_params.sigma_a_c = 8e-3;
  imu_params.sigma_bg = 0.03;
  imu_params.sigma_ba = 0.1;
  imu_params.sigma_gw_c = 4e-6;
  imu_params.sigma_aw_c = 4e-5;
  imu_params.tau = 3600;
  imu_params.g = 9.81007;
  imu_params.a0 = {0., 0., 0.};
  imu_params.rate = 200;

  pathgen::ImuMeasurementDeque imu_measurements;
  pathgen::ImuPathOptions path_options;
  path_options.add_noise_to_imu_measurements = false;
  std::vector<pathgen::PosePtr> imu_poses;
  std::vector<pathgen::PosePtr> poses;
  uint64_t duration = genPoses(&poses);
  path_options.duration = 1e-9*static_cast<double>(duration);

  pathgen::PoseSpline pose_spline =
    pathgen::PathGenerator::GenerateSplineFromPoses(poses, true);

  pathgen::IMUGenerator::GenerateInertialMeasurements(pose_spline, imu_params,
      path_options, &imu_measurements, &imu_poses);

  // ground truth acceleration and gyroscope calculation:
  // pose is traveling in circle of radius 1, through arc of size pi/3, in 2
  // seconds; rotational velocity is pi/6 rad/sec, centripetal acceleration is
  // v^2/r = (pi/6)^2, gravity is 9.81
  Eigen::Vector3d accel_gt(0, M_PI*M_PI/36., imu_params.g);
  Eigen::Vector3d gyro_gt(0., 0., -M_PI/6.);
  for (size_t t = 50; t < imu_measurements.size() - 50; t++) {
    if ((accel_gt - imu_measurements[t].measurement.accelerometers).norm() >
        1e-2) {
      std::cout << "Test failed: " << t << " " <<
        imu_measurements[t].measurement.accelerometers.transpose() << std::endl;
    }
    if ((gyro_gt - imu_measurements[t].measurement.gyroscopes).norm() > 1e-2) {
      std::cout << "Test failed: " << t << " " <<
        imu_measurements[t].measurement.gyroscopes.transpose() << std::endl;
    }
  }

  return 0;
}
