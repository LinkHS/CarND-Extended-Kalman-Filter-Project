#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  _x = VectorXd(4);

  // clang-format off
  // state covariance matrix
  _P = MatrixXd(4, 4);
  _P << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 100, 0,
        0, 0, 0, 100;

  // state transition matrix
  _F = MatrixXd(4, 4);
  _F << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  // measurement matrix - laser
  _H_laser_ = MatrixXd(2, 4);
  _H_laser_ << 1, 0, 0, 0,
               0, 1, 0, 0;

  //measurement covariance matrix - laser
  _R_laser_ = MatrixXd(2, 2);
  _R_laser_ << 0.0225, 0,
               0, 0.0225;

  // measurement matrix - radar
  _Hj_radar = MatrixXd(3, 4);

  //measurement covariance matrix - radar
  _R_radar_ = MatrixXd(3, 3);
  _R_radar_ << 0.09, 0, 0,
               0, 0.0009, 0,
               0, 0, 0.09;

  // process covariance matrix
  _Q = MatrixXd(4, 4);

  _I = MatrixXd::Identity(_x.size(), _x.size());

  // set the acceleration noise components
  _noise_ax = 9;
  _noise_ay = 9;

  // clang-format on
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::Predict(float delta_t)
{
  // std::cout << __func__ << "+" << std::endl;

  float dt_2 = delta_t * delta_t;
  float dt_3 = dt_2 * delta_t;
  float dt_4 = dt_3 * delta_t;

  // Update the state transition matrix F according to the new elapsed time.
  //   Time is measured in seconds.
  _F(0, 2) = delta_t;
  _F(1, 3) = delta_t;

  // Update the process noise covariance matrix.
  // clang-format off
  _Q << dt_4/4*_noise_ax, 0, dt_3/2*_noise_ax, 0,
        0, dt_4/4*_noise_ay, 0, dt_3/2*_noise_ay,
        dt_3/2*_noise_ax, 0, dt_2*_noise_ax, 0,
        0, dt_3/2*_noise_ay, 0, dt_2*_noise_ay;
  // clang-format on

  // std::cout << __func__ << "-> Predict" << std::endl;
  KalmanFilter::Predict();

  // std::cout << __func__ << "-" << std::endl;
}

void FusionEKF::Update(const MeasurementPackage &measurement_pack)
{
  // std::cout << __func__ << "+" << std::endl;
  VectorXd z;

  /*
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    z = VectorXd(3, 1);
    z << measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1],
        measurement_pack.raw_measurements_[2];

    // Recover state parameters
    VectorXd z_pred = VectorXd(3);
    ConvertvCartesianToPolar(_x, z_pred);
    VectorXd y = z - z_pred;

    // Normalize angle
    while (y(1) > M_PI)
      y(1) -= 2 * M_PI;
    while (y(1) < -M_PI)
      y(1) += 2 * M_PI;

    Calculate_Hj_radar(_x, &_Hj_radar);
    _p_R = &_R_radar_;
    _p_H = &_Hj_radar;

    // std::cout << __func__ << "-> Update kf" << std::endl;
    KalmanFilter::Update(nullptr, &y);
  }
  else
  {
    // Laser updates
    z = VectorXd(2, 1);
    z << measurement_pack.raw_measurements_[0],
        measurement_pack.raw_measurements_[1];

    _p_R = &_R_laser_;
    _p_H = &_H_laser_;

    // std::cout << __func__ << "-> Update ekf" << std::endl;
    KalmanFilter::Update(&z);
  }

  // std::cout << __func__ << "-" << std::endl;
}

/** 将笛卡尔坐标转化为
 * @ cartesian_in
 */
void FusionEKF::ConvertvCartesianToPolar(const Eigen::VectorXd &cartesian_in, Eigen::VectorXd &polar_out)
{
  float px = cartesian_in(0);
  float py = cartesian_in(1);
  float vx = cartesian_in(2);
  float vy = cartesian_in(3);

  float px2 = px * px;
  float py2 = py * py;
  float sumSQR = px2 + py2;
  assert(sumSQR != 0);

  float rho = sqrt(sumSQR);
  float theta = atan2(py, px);
  float rho_dot = (px * vx + py * vy) / rho;

  assert(polar_out.size() == 3);
  polar_out << rho, theta, rho_dot;
}

int FusionEKF::Calculate_Hj_radar(const VectorXd &x_state, MatrixXd *p_Hj)
{
  /**
   * Calculate a Jacobian here.
   */
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px * px + py * py;
  float c2 = sqrt(c1);
  float c3 = (c1 * c2);

  // check division by zero
  if (fabs(c1) < 0.0001)
  {
    std::cerr << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return -1;
  }

  // compute the Jacobian matrix
  // clang-format off
  *p_Hj << (px/c2), (py/c2), 0, 0,
           -(py/c1), (px/c1), 0, 0,
           py * (vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  // clang-format on

  return 0;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /**
   * Initialization
   */
  if (!is_initialized_)
  {
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      cout << "EKF: RADAR" << endl;

      // Convert radar from polar to cartesian coordinates
      //   and initialize state.
      float range = measurement_pack.raw_measurements_[0];
      float angle = measurement_pack.raw_measurements_[1];
      float range_rate = measurement_pack.raw_measurements_[2];

      float px = range * cos(angle);
      float py = range * sin(angle);
      float vx = range_rate * cos(angle);
      float vy = range_rate * sin(angle);

      _x << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      cout << "EKF: LASER" << endl;

      // Initialize state.
      _x << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1],
          0,
          0;

      _P(2, 2) = 1000;
      _P(3, 3) = 1000;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    std::cout << "End of Init" << std::endl;
    return;
  }

  /**
   * Prediction
   */
  float delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  Predict(delta_t);

  /**
   * Update
   */
  Update(measurement_pack);

  // print the output
  cout << "x_ = " << _x << endl;
  cout << "P_ = " << _P << endl;
}
