#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;


VectorXd calc_carterisan_to_polar(const VectorXd &x_state) {

  float px = x_state[0];
  float py = x_state[1];  
  float vx = x_state[2];
  float vy = x_state[3];

  //float rho, phi, rho_dot;
  float rho = sqrt( pow(px,2) + pow(py,2)  );
  float phi = atan2(py, px);
  float rho_dot;
  if (rho < 0.0001) 
    rho_dot = 0.0;
  else
    rho_dot = ( px * vx + py * vy ) / rho;

  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rho_dot;

  return z_pred;

}



// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
 	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_ * x_;

	VectorXd y = z - z_pred;

  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
 
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  const float DoublePI = 2 * M_PI;
  VectorXd z_pred = calc_carterisan_to_polar(x_);
  VectorXd y = z - z_pred;

  // if y[1] is over PI, then extract 2PI 
  while( y(1) > M_PI ) {
    y(1) -= DoublePI;
  }
  // if y[1] is under PI, then add 2PI
  while( y(1) < -M_PI ) {
    y(1) += DoublePI;
  }
  // the above adjustment is stablized to caculate RMSE. 
  // so that estimation close to grand truth is outcome. 

  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;


}

