#include <pathgen/Utils.h>

namespace pathgen{

void SavePoses(std::vector<PosePtr>& poses, const std::string filename){
    std::ofstream poses_file(filename, std::ios_base::trunc);

    for(std::vector<PosePtr>::iterator it = poses.begin();
        it != poses.end(); ++it)
    {
        poses_file << (*it)->translation().transpose().format(pathgen::kLongCsvFmt)
                   << ", " <<  (*it)->unit_quaternion().w() << "," << (*it)->unit_quaternion().x() << "," << (*it)->unit_quaternion().y() << "," << (*it)->unit_quaternion().z() << std::endl;
    }
}

void SaveSpline(const int num_measurements,
                const PoseSpline& spline,
                const double duration){

    // write all the spline poses to a file, at imu-rate
    std::ofstream spline_file("spline.txt", std::ios_base::trunc);
    const double step_size = (double)1.0/num_measurements;
    for(double ii = 0; ii <= 1.0; ii += step_size)
    {
        // sample spline uniformlly across points
        PointType pt = spline.position_spline(ii);
        PointType accel = spline.position_spline.derivatives(ii, 2).col(2)/(duration*duration);
        PointType vel = spline.position_spline.derivatives(ii, 2).col(1)/duration;

        Eigen::Quaterniond q;
        Eigen::Vector3d gyro;

        // use quaternion spline...perferred method
        q = spline.quat_spline(ii);
        gyro = (spline.quat_spline.derivative(ii)/duration);

        spline_file << ii*duration << ", "; // step size, scaled by the total duration
        spline_file << pt[0] << "," << pt[1] << "," << pt[2] << ","; //write point out
        spline_file << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << ", "; // quaternion
        spline_file << accel.transpose().format(pathgen::kLongCsvFmt) << ", " << // acceleration
                       vel.transpose().format(pathgen::kLongCsvFmt) << ", "; // velocity
        spline_file << gyro.transpose().format(pathgen::kLongCsvFmt) << std::endl; // angular velocity

    }

}

void SaveCameraPoses(const PoseSpline& pose_spline, const int camera_rate,
                     const int duration,
                     const double start_time){
    std::ofstream cam_poses_file("cam_poses.txt", std::ios_base::trunc);
    const size_t NUM_FRAMES = duration*camera_rate;

    for (size_t k = 0; k < NUM_FRAMES + 1; ++k) {
        double t = (double)1.0/NUM_FRAMES * (double)k;
        Pose T_WS = pose_spline(t);

        double frame_timestamp = start_time +
                (double(k) * (duration) / double(NUM_FRAMES));

        // write out pose to a file
        Eigen::Vector3d pqr = R2pqr(T_WS.rotationMatrix());
        cam_poses_file << frame_timestamp << ", " <<
                          T_WS.translation().transpose().format(pathgen::kLongCsvFmt) <<
                          ", " <<  pqr[0] << "," << pqr[1]<< "," << pqr[2] <<
                          std::endl;

    }
}
}
