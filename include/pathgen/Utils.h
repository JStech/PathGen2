#pragma once
#ifndef INCLUDE_PATHGEN_UTILS_H_
#define INCLUDE_PATHGEN_UTILS_H_
#include <Eigen/Eigen>

#include <pathgen/Types.h>
#include <pathgen/QuaternionSpline.h>
#include <pathgen/PoseSpline.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>

namespace pathgen {


static Eigen::IOFormat kCleanFmt(4, 0, ", ", ";\n", "", "");
static Eigen::IOFormat kLongFmt(Eigen::FullPrecision, 0, ", ", ";\n", "", "");
static Eigen::IOFormat kLongCsvFmt(Eigen::FullPrecision, 0, ", ", "\n", "", "");



inline Eigen::Matrix3d pqr2R(double p, double q, double r){
    double cx = cos(p), sx = sin(p);
    double cy = cos(q), sy = sin(q);
    double cz = cos(r), sz = sin(r);

    double ss_ = sx * sy;
    double cs_ = cx * sy;
    double sc_ = sx * cy;
    double cc_ = cx * cy;
    double c_s = cx * sz;
    double s_s = sx * sz;
    double cs__ = cy * sz;
    double cc__ = cy * cz;
    double s_c = sx * cz;
    double c_c = cx * cz;
    double ssc = ss_ * cz;
    double csc = cs_ * cz;
    double sss = ss_ * sz;
    double css = cs_ * sz;

    Eigen::Matrix3d R;
    R  <<  cc__, -c_s + ssc,  s_s + csc,
            cs__,  c_c + sss, -s_c + css,
            -sy,        sc_,        cc_;

    return R;
}

inline Eigen::Vector3d R2pqr(Eigen::Matrix3d R){

    Eigen::Vector3d euler;
    euler(0) = atan2( R(2,1), R(2,2) );   //roll, psi or p
    euler(1) = -asin( R(2,0) );           // pitch, theta or q
    euler(2) = atan2( R(1,0), R(0,0) );   // yaw, phi or r
    return euler;
}

inline Eigen::Matrix3d quat2R(Eigen::Quaterniond q){

    Eigen::Matrix3d R;
    R(0,0) = q.w()*q.w() + q.vec()[0]*q.vec()[0] - q.vec()[1]*q.vec()[1] - q.vec()[2]*q.vec()[2];
    R(0,1) = 2*( q.vec()[0]*q.vec()[1] - q.w()*q.vec()[2]);
    R(0,2) = 2*( q.vec()[0]*q.vec()[2] + q.w()*q.vec()[1]);

    R(1,0) = 2*( q.vec()[1]*q.vec()[0] + q.w()*q.vec()[2]);
    R(1,1) = q.w()*q.w() - q.vec()[0]*q.vec()[0] + q.vec()[1]*q.vec()[1] - q.vec()[2]*q.vec()[2];
    R(1,2) = 2*( q.vec()[1]*q.vec()[2] - q.w()*q.vec()[0]);

    R(2,0) = 2*( q.vec()[2]*q.vec()[0] - q.w()*q.vec()[1]);
    R(2,1) = 2*( q.vec()[2]*q.vec()[1] + q.w()*q.vec()[0]);
    R(2,2) = q.w()*q.w() - q.vec()[0]*q.vec()[0] - q.vec()[1]*q.vec()[1] + q.vec()[2]*q.vec()[2];

    return R;
}

inline double angleWrap(const double a)
{

    double angle = a;
    while( angle > M_PI ){
        angle = angle - 2*M_PI;
    }
    while( angle < -M_PI ){
        angle = angle + 2*M_PI;
    }

    return angle;
}

inline void angleWrapPQR(Eigen::Vector3d* pqr)
{
    (*pqr)(0) = angleWrap( (*pqr)(0) );
    (*pqr)(1) = angleWrap( (*pqr)(1) );
    (*pqr)(2) = angleWrap( (*pqr)(2) );
}


inline Eigen::Quaterniond euler2quat(Eigen::Vector3d euler) // euler angles in roll, pitch, yaw order
{
    angleWrapPQR(&euler);

    double cosX = cos(euler(0)/2.0);
    double cosY = cos(euler(1)/2.0);
    double cosZ = cos(euler(2)/2.0);

    double sinX = sin(euler(0)/2.0);
    double sinY = sin(euler(1)/2.0);
    double sinZ = sin(euler(2)/2.0);

    Eigen::Vector4d coeffs(cosX*cosY*cosZ + sinX*sinY*sinZ,
                           sinX*cosY*cosZ - cosX*sinY*sinZ,
                           cosX*sinY*cosZ + sinX*cosY*sinZ,
                           cosX*cosY*sinZ - sinX*sinY*cosZ);

    return Eigen::Quaterniond(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
}

void SavePoses(std::vector<PosePtr>& poses, const std::string filename);
void SaveCameraPoses(const PoseSpline& pose_spline, const int camera_rate,
                     const int duration,
                     const double start_time = 0);
void SaveSpline(const int num_measurements,
                const PoseSpline& spline,
                const double duration);



}

#endif
