#include "library.h"

double Fusion::predict(const VectorXd &u, const MatrixXd &Q) {
    //Assuming x is divided into position and velocity
    assert(N % 2 == 0);

    //Check input size
    assert(u.rows() == N && Q.rows() == N && Q.cols() == N);

    //Compute dt
    double dt = getDt();

    //Compute F matrix for constant velocity movement
    MatrixXd F(N, N);
    F << I(N / 2), I(N / 2) * dt,
            Z(N / 2), I(N / 2);

    //Predict according to kalman equations
    x = F * x + u*dt;
    P = F * P * F.transpose() + Q*dt;

    return dt;
}


void Fusion::update(const VectorXd &z, const MatrixXd &R) {
    //Check input size
    assert(z.rows() == N / 2 && R.cols() == N / 2 && R.rows() == N / 2);

    //Compute H matrix so that we only observe position, not velocity
    MatrixXd H(N / 2, N);
    H << I(N / 2), Z(N / 2);

    //Update according to kalman equations
    VectorXd y = z - H * x;
    MatrixXd S = H * P * H.transpose() + R;
    MatrixXd K = P * H.transpose() * S.inverse();
    x = x + K * y;
    P = (I(N) - K * H) * P;
}
