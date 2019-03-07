#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const MatrixXd &x,
                        const MatrixXd &P,
                        const MatrixXd &F,
                        const MatrixXd &Q,
                        const MatrixXd *p_H,
                        const MatrixXd *p_R)
{
  _x = x;
  _P = P;
  _F = F;
  _Q = Q;
  // clang-format off
  if(p_H != nullptr) {_p_H = p_H;};
  if(p_R != nullptr) {_p_R = p_R;};
  // clang-format on

  _I = MatrixXd::Identity(_x.size(), _x.size());
}

void KalmanFilter::Predict(void)
{
  _x = _F * _x;
  _P = _F * _P * _F.transpose() + _Q;
}

/**
 * Update the state by using Kalman Filter equations
 * 
 * @p_z: For extended kalman filter. 
 *           Set to 'nullptr' for normal kalman filter.
 */
void KalmanFilter::Update(const VectorXd *p_z, const VectorXd *p_y)
{
  VectorXd y;

  if (p_z != nullptr)
  {
    // std::cout << __func__ << "1" << std::endl;
    VectorXd z_pred = *_p_H * _x;
    // std::cout << __func__ << "2" << std::endl;
    y = *p_z - z_pred;
  }
  else if (p_y != nullptr) // For extended kalman filter
  {
    // std::cout << __func__ << "12" << std::endl;
    y = *p_y;
  }
  else
  {
    std::cerr << "Unsupported arguments" << std::endl;
  }

  // std::cout << __func__ << "3" << std::endl;
  MatrixXd Ht = (*_p_H).transpose();

  // std::cout << __func__ << "4" << std::endl;
  MatrixXd S = *_p_H * _P * Ht + *_p_R;

  // std::cout << __func__ << "5" << std::endl;
  MatrixXd Si = S.inverse();

  // std::cout << __func__ << "6" << std::endl;
  MatrixXd PHt = _P * Ht;

  // std::cout << __func__ << "7" << std::endl;
  MatrixXd K = PHt * Si;

  // std::cout << __func__ << "8" << std::endl;
  _x = _x + (K * y);

  // std::cout << __func__ << "9" << std::endl;
  _P = (_I - K * (*_p_H)) * _P;

  // std::cout << __func__ << "-" << std::endl;
}