#include<iostream>
#include "library.h"
#include "multivar_noise.h"
#include <fstream>
#include <cmath>
#include <cstdlib> // for strtol()
#include <utility>

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

VectorXd vehicle_state_transition_function(const VectorXd &x, const VectorXd &u) {
    double m = 1292.2;
    double a = 1.006;
    double b = 1.534;
    double I = 2380.7;

    //Construct F
    VectorXd f(6);
    f << x(5) * cos(x(2)) - x(4) * sin(x(2)),
            x(4) * cos(x(2)) + x(5) * sin(x(2)),
            x(3),
            (b * u(0) * u(1) * 2) / (I * (u(1) + 1)),
            -x(3) * x(5),
            -x(3) * x(5) + (u(0) - (u(0) * pow(u(1), 2)) / (u(1) + 1)) / m;

    return f;
};

MatrixXd vehicle_state_transition_jacobian(const VectorXd &x, const VectorXd &u) {
    //Construct J
    MatrixXd J = MatrixXd::Zero(6, 6);
    J(0, 2) = -x(4) * cos(x(2)) - x(5) * sin(x(2));
    J(0, 4) = -sin(x(2));
    J(0, 5) = cos(x(2));
    J(1, 2) = x(5) * cos(x(2)) - x(4) * sin(x(2));
    J(1, 4) = cos(x(2));
    J(1, 5) = sin(x(2));
    J(2, 3) = 1;
    J(4, 3) = -x(5);
    J(4, 5) = -x(3);
    J(5, 3) = -x(5);
    J(5, 5) = -x(3);

    return J;
};

typedef VectorXd state_type;
runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper;
struct ode {
    Fusion *f;
    VectorXd u;
    normal_random_variable v;

    ode(Fusion *f, VectorXd u, normal_random_variable v) : f(f), u(move(u)), v(move(v)) {}

    void operator()(state_type const &x, state_type &dxdt, double t) const {
        dxdt = f->state_transition_function(x, u)+v();
    }
};


int main(int argc, char *argv[]) {
    int testCaseIndex = 3;
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
            dt = 0.03;
            T = 4;
            f = Fusion(N);
            f.setConstantDt(dt);
            f.state_transition_function = vehicle_state_transition_function;
            f.state_transition_jacobian = vehicle_state_transition_jacobian;
            x = VectorXd::Zero(N);
            Q = MatrixXd::Random(N, N) / 5;
            R = MatrixXd::Random(N / 2, N / 2) / 10;
            u = VectorXd(2);
            u << 5000, 5 * 3.14 / 180.0; //Pr, d
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
        if (testCaseIndex == 3 && i * dt >= T / 4) {
            u(1) = -abs(u(1));
        }
        integrate_const(stepper, ode(&f, u, v), x, 0.0, dt, dt / 100);
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
