#include<iostream>
#include "library/library.h"
#include "multivar_noise.h"
#include <fstream>
#include <cmath>
#include <cstdlib> // for strtol()
#include <utility>
#include "functions/functions.h"

using namespace std;


int main(int argc, char *argv[]) {
    //Input variables
    int testCaseIndex = 3;
    if (argc >= 2) {
        testCaseIndex = strtol(argv[1], nullptr, 10);
    }
    bool only_ground_truth = false;
    if (argc >= 3) {
        only_ground_truth = strtol(argv[2], nullptr, 10) == 1;
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
    MatrixXd R2 = MatrixXd::Identity(2, 2);
    int use_R2_every_x_steps = 0;
    double max_steering_angle = 0;
    switch (testCaseIndex) {
        case 1: //constant_movement
            N = 4;
            dt = 0.03;
            T = 1.0;
            f = Fusion(N);
            f.setConstantDt(dt);
            x = 10 * VectorXd::Random(N);
            Q = MatrixXd::Random(N, N) / 10;
            R = MatrixXd::Random(N / 2, N / 2) / 10;
            u = VectorXd(N);
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
            x = VectorXd(N);
            x << 1.033366313746765, 0, 0, -.05849376854515592;
            Q = MatrixXd::Random(N, N) / 1000;
            R = MatrixXd::Random(N / 2, N / 2) / 1000;
            u = VectorXd(N);
            u << 0, 0, 0, 0;
            break;
        case 3: //Vehicle dynamics
            N = 6;
            dt = 0.3;
            T = 100;
            x = VectorXd::Zero(N);
            x(3) = 1; // setting vn to some value (if vn=zero the slip angle is nan)
            // Also setting the P to zero to let it know it is a perfect guess
            // (if not it can go to nan when guessing the answer at the begining)
            u = VectorXd(N);
            u << 0, 0, 0, 0, 0, 0;
            f = Fusion(N, x, DiagonalMatrix<double, 6>(u));
            f.setConstantDt(dt);
            f.state_transition_function = vehicle_state_transition_function;
            f.state_transition_jacobian = vehicle_state_transition_jacobian;
            f.observation_function = vehicle_observation_function;
            f.observation_jacobian = vehicle_observation_jacobian;
            Q = MatrixXd::Identity(6, 6) / 20;
            R = DiagonalMatrix<double, 3>(1. / 20, 1. / 20, 1. / 20);
            R2 = DiagonalMatrix<double, 3>(1. / 20, 1. / 20, 1. / 20);
            use_R2_every_x_steps = 0;
            u = VectorXd(3);
            max_steering_angle = 10 * EIGEN_PI / 180.0;
            u << 0, 0, 0; //Pf, Pr, d

            break;
        default:
            exit(1);
    }



    //Make Q and R into positive semi-definite matrices
    Q = Q.transpose() * Q;
    R = R.transpose() * R;
    MatrixXd Rf = MatrixXd::Zero(N, N);
    Rf.bottomRightCorner(N / 2, N / 2) = R;
    R2 = R2.transpose() * R2;
    MatrixXd R2f = MatrixXd::Zero(N, N);
    R2f.bottomRightCorner(N / 2, N / 2) = R;

    //Initialize noise generators for v and w
    normal_random_variable v{Q};
    normal_random_variable w{R};
    normal_random_variable w2{R2};
    normal_random_variable zero_noise{MatrixXd::Zero(N, N)};

    //Start simulation
    VectorXd ground_truth = x;
    for (int i = 0; i < (int) T / dt; i++) {
        //Simulate movement with simulated process noise and try to predict it
        if (testCaseIndex == 3 and i % 10 == 0) {
            u(2) += fRand(-1., 1.) * 5 * EIGEN_PI / 180; // NOLINT(cert-msc30-c,cert-msc50-cpp)
            u(2) = clamp(u(2), -max_steering_angle, max_steering_angle);
        }
        VectorXd truth_prev = ground_truth;
        VectorXd x_prev = x;
        integrate(dt, f.state_transition_function, ground_truth, u, v);
        integrate(dt, f.state_transition_function, x, u, v);

        VectorXd z = f.observation_function(ground_truth);
        if (testCaseIndex == 3) { //should be a boolean like use_deltas above
            VectorXd vars_dot = (ground_truth.head(N / 2) - truth_prev.head(N / 2)) / dt;
            double th = truth_prev(2);
            z << vars_dot(0) * cos(th) + vars_dot(1) * sin(th), //Vn
                    -vars_dot(0) * sin(th) + vars_dot(1) * cos(th), //Ve
                    vars_dot(2); //th_dot
        }

        bool use_R2 = use_R2_every_x_steps > 0 && i % use_R2_every_x_steps == 0 && i > 0;
        if (use_R2)
            z += w2();
        else
            z += w();

        if (testCaseIndex == 3) {//should be a boolean like use_deltas above
            double th = x_prev(2);
            VectorXd vars_dot(3);
            vars_dot << z(0) * cos(th) - z(1) * sin(th),
                    z(0) * sin(th) + z(1) * cos(th), z(2);
            x.head(N / 2) = vars_dot * dt + x_prev.head(N / 2);
        }
        if (!only_ground_truth) {
            if (isnan(f.getP()(1, 1))) cout << "Is nan before predict" << endl;
            f.predict(u, Q);
            if (isnan(f.getP()(1, 1))) cout << "Is nan after predict" << endl;

            //Update kalman filter with simulated measurement noise
            if (use_R2)
                f.update(z, R2);
            else
                f.update(z, R);
            if (isnan(f.getP()(1, 1))) cout << "Is nan after update" << endl;
        }
        //Print information
        if (isnan(f.getP()(1, 1))) {
            cout << "isnan i=" << i << endl;
            return 0;
        }
        file << ground_truth.head(2).format(csv) << "," << x.head(2).format(csv) << ","
             << f.getx().head(2).format(csv)
             << "," << f.getP().topLeftCorner(2, 2).format(csv) << endl;
    }
    file.close();
}
