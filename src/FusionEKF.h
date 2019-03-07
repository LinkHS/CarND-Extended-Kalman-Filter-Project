#ifndef FusionEKF_H_
#define FusionEKF_H_

#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "kalman_filter.h"
#include "measurement_package.h"
#include "tools.h"

class FusionEKF : public KalmanFilter
{
public:
  /**
   * Constructor.
   */
  FusionEKF();

  /**
   * Destructor.
   */
  virtual ~FusionEKF();

  /**
   * Run the whole flow of the Kalman Filter from here.
   */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  void Predict(float delta_t);

  void Update(const MeasurementPackage &measurement_pack);

  int Calculate_Hj_radar(const Eigen::VectorXd &x_state, Eigen::MatrixXd *p_Hj);

  void ConvertvCartesianToPolar(const Eigen::VectorXd &cartesian_in, Eigen::VectorXd &polar_out);

  void GetState(double *p_px, double *p_py, double *p_vx, double *p_vy)
  {
    *p_px = _x(0);
    *p_py = _x(1);
    *p_vx = _x(2);
    *p_vy = _x(3);
  }

private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  Eigen::MatrixXd _R_laser_;
  Eigen::MatrixXd _R_radar_;
  Eigen::MatrixXd _H_laser_;
  Eigen::MatrixXd _Hj_radar;

  float _noise_ax;
  float _noise_ay;
};

#endif // FusionEKF_H_
