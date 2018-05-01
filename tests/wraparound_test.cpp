#include <pathgen/PathGen.h>
#include <Eigen/Geometry>
#include <glog/logging.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define TIME 500000000
#define STEP  10000000

uint64_t genPoses(std::vector<pathgen::PosePtr>* poses) {
  Eigen::Matrix4d T;
  Eigen::Quaterniond q(1, 0, 0, 0);

  for (uint64_t ts = 0; ts < TIME; ts += STEP) {
    double theta = static_cast<double>(ts)/static_cast<double>(TIME) * M_PI;
    T.topRightCorner(3, 1) = Eigen::Vector3d(
        sin(M_PI/2 + theta), cos(M_PI/2 + theta), 0);
    q.w() = cos((M_PI/2 + theta)/2);
    q.x() = 0;
    q.y() = 0;
    q.z() = sin((M_PI/2 + theta)/2);
    T.topLeftCorner(3, 3) = q.toRotationMatrix();
    pathgen::PosePtr new_pose(new pathgen::Pose(T));
    poses->push_back(new_pose);
  }
  return TIME;
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
  path_options.add_noise_to_imu_measurements = true;
  std::vector<pathgen::PosePtr> imu_poses;
  std::vector<pathgen::PosePtr> poses;
  uint64_t duration = genPoses(&poses);
  path_options.duration = 1e-9*static_cast<double>(duration);

  pathgen::PoseSpline pose_spline =
    pathgen::PathGenerator::GenerateSplineFromPoses(poses, true);

  pathgen::IMUGenerator::GenerateInertialMeasurements(pose_spline, imu_params,
      path_options, &imu_measurements, &imu_poses);

  // TODO test results somehow, rather than just printing them
  pathgen::IMUGenerator::SaveMeasurementFiles(imu_measurements);
  pathgen::SavePoses(imu_poses, "imu_poses.txt");

  return 0;
}
