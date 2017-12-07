#pragma once
#ifndef PATHGEN_IMU_GENERATOR_H
#define PATHGEN_IMU_GENERATOR_H

#include <pathgen/Types.h>
#include <pathgen/Measurements.h>
#include <pathgen/PoseSpline.h>
#include <pathgen/Utils.h>
#include <Eigen/Core>
#include <iostream>
#include <fstream>

namespace pathgen{

struct ImuPathOptions{
    ImuPathOptions(){}

    double duration = 10.0;
    bool add_noise_to_imu_measurements = true;

};

class IMUGenerator
{
public:
    IMUGenerator(){}

    static void GenerateInertialMeasurements(const PoseSpline& pose_spline,
            const ImuParameters& imu_parameters,
            const ImuPathOptions &path_options,
            ImuMeasurementDeque *imu_measurements,
            std::vector<PosePtr> *imu_poses);

    static void SaveMeasurementFiles(const ImuMeasurementDeque& meas);

};

}
#endif
