//
// Created by Inigo on 23/08/2020.
//

#include "functions.h"

#include <utility>

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
    double m = 1292.2;       //Vehicle mass
    double I = 2380.7;       //Vehicle inertia
    double a = 1.006;        //CG to front axle
    double b = 1.534;        //CG to rear axle
    double c = 20000;        //cornering stiffness

    //Construct F
    VectorXd f(6);
    f << x(5) * cos(x(2)) - x(4) * sin(x(2)),
            x(4) * cos(x(2)) + x(5) * sin(x(2)),
            x(3),
            (a * c * (u(2) - (x(4) + a * x(3)) / x(5)) + (b * c * (x(4) - b * x(3))) / x(5)) / I,
            -x(3) * x(5) - ((c * (x(4) - b * x(3))) / x(5) - c * cos(u(2)) * (u(2) - (x(4) + a * x(3)) / x(5))) / m,
            x(3) * x(4) + (u(1) - c * sin(u(2)) * (u(2) - (x(4) + a * x(3)) / x(5))) / m;


    return f;
};

MatrixXd vehicle_state_transition_jacobian(const VectorXd &x, const VectorXd &u) {
    double m = 1292.2;       //Vehicle mass
    double I = 2380.7;       //Vehicle inertia
    double a = 1.006 / 2;        //CG to front axle
    double b = 1.534 / 2;        //CG to rear axle
    double c = 20000;        //cornering stiffness
    //Construct J
    MatrixXd J = MatrixXd::Zero(6, 6);

    J(0, 2) = -x(4) * cos(x(2)) - x(5) * sin(x(2));
    J(0, 4) = -sin(x(2));
    J(0, 5) = cos(x(2));
    J(1, 2) = x(5) * cos(x(2)) - x(4) * sin(x(2));
    J(1, 4) = cos(x(2));
    J(1, 5) = sin(x(2));
    J(2, 3) = 1.0;
    J(3, 3) = -(((a * a) * c) / x(5) + ((b * b) * c) / x(5)) / I;
    J(3, 4) = -((a * c) / x(5) - (b * c) / x(5)) / I;
    J(3, 5) = (a * c * 1.0 / pow(x(5), 2) * (x(4) + a * x(3)) - b * c * 1.0 / pow(x(5), 2) * (x(4) - b * x(3))) / I;
    J(4, 3) = -x(5) + ((b * c) / x(5) - (a * c * cos(u(2))) / x(5)) / m;
    J(4, 4) = -(c / x(5) + (c * cos(u(2))) / x(5)) / m;
    J(4, 5) = -x(3) +
              (c * 1.0 / pow(x(5), 2) * (x(4) - b * x(3)) + c * 1.0 / pow(x(5), 2) * cos(u(2)) * (x(4) + a * x(3))) / m;
    J(5, 3) = x(4) + (a * c * sin(u(2))) / (m * x(5));
    J(5, 4) = x(3) + (c * sin(u(2))) / (m * x(5));
    J(5, 5) = -(c * 1.0 / pow(x(5), 2) * sin(u(2)) * (x(4) + a * x(3))) / m;
    return J;
};


runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper2;

void
integrate(double dt, const function<VectorXd(VectorXd const &, VectorXd const &)> &f, VectorXd &x, const VectorXd &u,
          const normal_random_variable &v) {
    integrate_const(stepper2, ode(f, u, v), x, 0.0, dt, dt / 100);
}

double fRand(double fMin, double fMax) {
    double f = (double) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}