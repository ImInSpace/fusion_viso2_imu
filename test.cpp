#include<iostream>
#include "library.h"
#include "multivar_noise.h"
#include <fstream>
#include <cmath>
#include <cstdlib> // for strtol()

using namespace std;


MatrixXd RTBP_state_transition_jacobian(const VectorXd &x, const VectorXd &u) {
    double mu = pow(10, -2);
    //Compute two radii
    double r1 = sqrt(pow(x(1), 2) + pow(x(0) - mu, 2));
    double r2 = sqrt(pow(x(1), 2) + pow(x(0) - mu + 1, 2));

    //Construct Jacobian Matrix
    MatrixXd J(4, 4);
    J.row(0) << 0, 0, 1, 0;
    J.row(1) << 0, 0, 0, 1;
    J(2, 0) = (mu - 1) / pow(r1, 3) - mu / pow(r2, 3) - (3 * (mu - 1) * pow(2 * x(0) - 2 * mu, 2)) / (4 * pow(r1, 5)) +
              (3 * mu * pow(2 * x(0) - 2 * mu + 2, 2)) / (4 * pow(r2, 5)) + 1;
    J(2, 1) = J(3, 0) = (3 * x(1) * mu * (2 * x(0) - 2 * mu + 2)) / (2 * pow(r2, 5)) -
                        (3 * x(1) * (mu - 1) * (2 * x(0) - 2 * mu)) / (2 * pow(r1, 5));
    J(3, 1) = (mu - 1) / pow(r1, 3) - mu / pow(r2, 3) + (3 * pow(x(1), 2) * mu) / pow(r2, 5) -
              (3 * pow(x(1), 2) * (mu - 1)) / pow(r1, 5) + 1;
    J.bottomRightCorner(2, 2) << 0, 2, -2, 0;

    return J;
};

VectorXd RTBP_state_transition_function(const VectorXd &x, const VectorXd &u) {
    double mu = pow(10, -2);
    //Compute two radii
    double r1 = sqrt(pow(x(1), 2) + pow(x(0) - mu, 2));
    double r2 = sqrt(pow(x(1), 2) + pow(x(0) - mu + 1, 2));

    //Construct F
    VectorXd f(4);
    f.head(2) = x.tail(2);
    f(2) = x(0) + 2 * x(3) - (mu * (2 * x(0) - 2 * mu + 2)) / (2 * pow(r2, 3)) +
           ((mu - 1) * (2 * x(0) - 2 * mu)) / (2 * pow(r1, 3));
    f(3) = x(1) - 2 * x(2) + (x(1) * (mu - 1)) / pow(r1, 3) - (x(1) * mu) / pow(r2, 3);

    return f;
};

int main(int argc, char *argv[]) {
    int testCaseIndex = 1;
    if (argc >= 2) {
        testCaseIndex = strtol(argv[1], nullptr, 10);
    }
    //Setup output format and file
    IOFormat singleLine(StreamPrecision, DontAlignCols, ",\t", ";\t", "", "", "[", "]");
    IOFormat csv(FullPrecision, DontAlignCols, ",", ",", "", "", "", "");
    ofstream file;
    file.open("data.csv", ios::trunc);

    //Setup random
    srand(time(NULL));


    int N;
    double dt;
    double T;
    Fusion f(0);
    VectorXd x;
    VectorXd u;
    MatrixXd Q;
    MatrixXd R;
    switch (testCaseIndex) {
        case 1: //constant_movement
            N = 4;
            dt = 0.03;
            T = 1.0;
            f = Fusion(N);
            f.setConstantDt(dt);
            x = 10 * VectorXd::Random(4);
            Q = MatrixXd::Random(4, 4) / 10;
            R = MatrixXd::Random(2, 2) / 10;
            u=VectorXd(4);
            u << 0, 0, 0, -0.1;
            break;
        case 2: //RTBP
            N = 4;
            dt = 0.03;
            T = 6.2296051135204102;
            f = Fusion(N);
            f.setConstantDt(dt);
            f.state_transition_function = RTBP_state_transition_function;
            f.state_transition_jacobian = RTBP_state_transition_jacobian;

            x = VectorXd(4);
            x << 1.033366313746765, 0, 0, -.05849376854515592;
            Q = MatrixXd::Random(4, 4) / 1000;
            R = MatrixXd::Random(2, 2) / 1000;
            u=VectorXd(4);
            u << 0, 0, 0, 0;
            break;
        default:
            exit(1);
    }


    //Make Q and R into positive semi-definite matrices
    Q = Q.transpose() * Q;
    R = R.transpose() * R;

    //Initialize noise generators for v and w
    normal_random_variable v{Q};
    normal_random_variable w{R};

    //Start simulation
    bool isnanprinted = false;
    for (int i = 0; i < (int) T / dt; i++) {
        //Simulate movement with simulated process noise and try to predict it
        x += f.state_transition_function(x, u) * dt;
        f.predict(u, Q);

        //Update kalman filter with simulated measurement noise
        VectorXd z = f.observation_function(x) + w();
        f.update(z, R);

        //Print information
        if (isnan(f.getP()(1, 1))) {
            if (!isnanprinted) {
                cout << "isnan i=" << i << endl;
                isnanprinted = true;
            }
            file << x.head(2).format(csv) << "," << z.head(2).format(csv) << ","
                 << "-1,-1"
                 << "," << "-1,-1,-1,-1" << endl;
        } else {
            file << x.head(2).format(csv) << "," << z.head(2).format(csv) << ","
                 << f.getx().head(2).format(csv)
                 << "," << f.getP().topLeftCorner(2, 2).format(csv) << endl;
        }
    }
    file.close();
}
