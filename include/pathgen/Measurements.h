#pragma once
#ifndef INCLUDE_PATHGEN_MEASUREMENTS_H
#define INCLUDE_PATHGEN_MEASUREMENTS_H

#include <deque>
#include <vector>
#include <memory>
#include <Eigen/Dense>

namespace pathgen {

struct ImuParameters{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    ImuParameters(){
        // set default parameters
        a0.setZero();
        g = 9.81;
        a_max = 1000.0;
        g_max = 1000.0;
        rate = 100;  // Hz
        sigma_g_c = 6.0e-4;
        sigma_a_c = 2.0e-3;
        sigma_gw_c = 3.0e-6;
        sigma_aw_c = 2.0e-5;
        sigma_bg = 0.03;
        sigma_ba = 0.1;
        tau = 3600.0;
    }

    double a_max;  ///< Accelerometer saturation. [m/s^2]
    double g_max;  ///< Gyroscope saturation. [rad/s]
    double sigma_g_c;  ///< Gyroscope noise density. [rad/s/sqrt(Hz)]
    double sigma_bg;  ///< Initial gyroscope bias. [rad/s]
    double sigma_a_c;  ///< Accelerometer noise density. [m/s^2/sqrt(Hz)]
    double sigma_ba;  ///< Initial accelerometer bias [m/s^2]
    double sigma_gw_c; ///< Gyroscope drift noise density. [rad/s^s/sqrt(Hz)]
    double sigma_aw_c; ///< Accelerometer drift noise density. [m/s^2/sqrt(Hz)]
    double tau;  ///< Reversion time constant of accerometer bias. [s]
    double g;  ///< Earth acceleration.
    Eigen::Vector3d a0;  ///< Mean of the prior accelerometer bias.
    int rate;  ///< IMU rate in Hz.
};

/**
 * \brief Generic measurements
 *
 * They always come with a timestamp such that we can perform
 * any kind of asynchronous operation.
 * \tparam MEASUREMENT_T Measurement data type.
 */
template<class MEASUREMENT_T>
struct Measurement {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double timestamp;      ///< Measurement timestamp
    MEASUREMENT_T measurement;  ///< Actual measurement.
    int sensor_id = -1;          ///< Sensor ID. E.g. camera index in a multicamera setup

    /// \brief Default constructor.
    Measurement()
        : timestamp(0.0) {
    }
    /**
   * @brief Constructor
   * @param timeStamp_ Measurement timestamp.
   * @param measurement_ Actual measurement.
   * @param sensorId Sensor ID (optional).
   */
    Measurement(const double& timestamp_, const MEASUREMENT_T& measurement_,
                int sensor_id_ = -1)
        : timestamp(timestamp_),
          measurement(measurement_),
          sensor_id(sensor_id_) {
    }
};

/// \brief IMU measurements.
struct ImuSensorReadings {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// \brief Default constructor.
    ImuSensorReadings()
        : gyroscopes(),
          accelerometers() {
    }
    /**
   * @brief Constructor.
   * @param gyroscopes_ Gyroscope measurement.
   * @param accelerometers_ Accelerometer measurement.
   */
    ImuSensorReadings(Eigen::Vector3d gyroscopes_,
                      Eigen::Vector3d accelerometers_)
        : gyroscopes(gyroscopes_),
          accelerometers(accelerometers_) {
    }
    Eigen::Vector3d gyroscopes;     ///< Gyroscope measurement.
    Eigen::Vector3d accelerometers; ///< Accelerometer measurement.
};

typedef Measurement<ImuSensorReadings> ImuMeasurement;
typedef std::deque<ImuMeasurement, Eigen::aligned_allocator<ImuMeasurement> > ImuMeasurementDeque;
typedef std::deque<const ImuMeasurement, Eigen::aligned_allocator<ImuMeasurement> > ConstImuMeasurementDeque;

typedef Eigen::Matrix<double, 9, 1> SpeedAndBias;

}  // namespace pathgen

#endif // INCLUDE_PATHGEN_MEASUREMENTS_H
