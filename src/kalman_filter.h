#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"

class KalmanFilter
{
public:
  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param _x Initial state
   * @param _P Initial state covariance
   * @param _F Transition matrix
   * @param _H Measurement matrix
   * @param _R Measurement covariance matrix
   * @param _Q Process covariance matrix
   */
  virtual void Init(const Eigen::MatrixXd &x,
                    const Eigen::MatrixXd &P,
                    const Eigen::MatrixXd &F,
                    const Eigen::MatrixXd &Q,
                    const Eigen::MatrixXd *p_H = nullptr,
                    const Eigen::MatrixXd *p_R = nullptr);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  virtual void Predict(void);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  virtual void Update(const Eigen::VectorXd *p_z, const Eigen::VectorXd *p_zpred=nullptr);

protected:
  // state vector
  Eigen::VectorXd _x;

  // state covariance matrix
  Eigen::MatrixXd _P;

  // state transition matrix
  Eigen::MatrixXd _F;

  // process covariance matrix
  Eigen::MatrixXd _Q;

  // measurement matrix
  const Eigen::MatrixXd *_p_H;

  // measurement covariance matrix
  const Eigen::MatrixXd *_p_R;

  // identity matrix
  Eigen::MatrixXd _I;
};

#endif // KALMAN_FILTER_H_
