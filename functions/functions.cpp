//
// Created by Inigo on 23/08/2020.
//

#include "functions.h"

#include <utility>

MatrixXd RTBP_state_transition_jacobian(const VectorXd& x, const VectorXd& u)
{
    double mu = pow(10, -2);
    // Compute two radii
    double r1 = sqrt(pow(x(1), 2) + pow(x(0) - mu, 2));
    double r2 = sqrt(pow(x(1), 2) + pow(x(0) - mu + 1, 2));

    // Construct Jacobian Matrix
    MatrixXd J(4, 4);
    J.row(0) << 0, 0, 1, 0;
    J.row(1) << 0, 0, 0, 1;
    J(2, 0) = (mu - 1) / pow(r1, 3) - mu / pow(r2, 3)
              - (3 * (mu - 1) * pow(2 * x(0) - 2 * mu, 2)) / (4 * pow(r1, 5))
              + (3 * mu * pow(2 * x(0) - 2 * mu + 2, 2)) / (4 * pow(r2, 5)) + 1;
    J(2, 1) = J(3, 0) = (3 * x(1) * mu * (2 * x(0) - 2 * mu + 2)) / (2 * pow(r2, 5))
                        - (3 * x(1) * (mu - 1) * (2 * x(0) - 2 * mu)) / (2 * pow(r1, 5));
    J(3, 1) = (mu - 1) / pow(r1, 3) - mu / pow(r2, 3) + (3 * pow(x(1), 2) * mu) / pow(r2, 5)
              - (3 * pow(x(1), 2) * (mu - 1)) / pow(r1, 5) + 1;
    J.bottomRightCorner(2, 2) << 0, 2, -2, 0;

    return J;
}

VectorXd RTBP_state_transition_function(const VectorXd& x, const VectorXd& u)
{
    double mu = pow(10, -2);
    // Compute two radii
    double r1 = sqrt(pow(x(1), 2) + pow(x(0) - mu, 2));
    double r2 = sqrt(pow(x(1), 2) + pow(x(0) - mu + 1, 2));

    // Construct F
    VectorXd f(4);
    f.head(2) = x.tail(2);
    f(2) = x(0) + 2 * x(3) - (mu * (2 * x(0) - 2 * mu + 2)) / (2 * pow(r2, 3))
           + ((mu - 1) * (2 * x(0) - 2 * mu)) / (2 * pow(r1, 3));
    f(3) = x(1) - 2 * x(2) + (x(1) * (mu - 1)) / pow(r1, 3) - (x(1) * mu) / pow(r2, 3);

    return f;
}

VectorXd vehicle_state_transition_function(const VectorXd& x, const VectorXd& u)
{
    double m = 1292.2;  // Vehicle mass
    double I = 2380.7;  // Vehicle inertia
    double a = 1.006;   // CG to front axle
    double b = 1.534;   // CG to rear axle
    double c = 20000;   // cornering stiffness

    // Construct F
    VectorXd f(6);
    f << x(3) * cos(x(2)) - x(4) * sin(x(2)), x(4) * cos(x(2)) + x(3) * sin(x(2)), x(5),
        x(4) * x(5) + (u(1) - c * sin(u(2)) * (u(2) - (x(4) + a * x(5)) / x(3))) / m,
        -x(3) * x(5)
            - ((c * (x(4) - b * x(5))) / x(3) - c * cos(u(2)) * (u(2) - (x(4) + a * x(5)) / x(3)))
                  / m,
        (a * c * (u(2) - (x(4) + a * x(5)) / x(3)) + (b * c * (x(4) - b * x(5))) / x(3)) / I;

    return f;
}

MatrixXd vehicle_state_transition_jacobian(const VectorXd& x, const VectorXd& u)
{
    double m = 1292.2;     // Vehicle mass
    double I = 2380.7;     // Vehicle inertia
    double a = 1.006 / 2;  // CG to front axle
    double b = 1.534 / 2;  // CG to rear axle
    double c = 20000;      // cornering stiffness
    // Construct J
    MatrixXd J = MatrixXd::Zero(6, 6);

    J(0, 2) = -x(4) * cos(x(2)) - x(3) * sin(x(2));
    J(0, 3) = cos(x(2));
    J(0, 4) = -sin(x(2));
    J(1, 2) = x(3) * cos(x(2)) - x(4) * sin(x(2));
    J(1, 3) = sin(x(2));
    J(1, 4) = cos(x(2));
    J(2, 5) = 1.0;
    J(3, 3) = -(c * 1.0 / pow(x(3), 2) * sin(u(2)) * (x(4) + a * x(5))) / m;
    J(3, 4) = x(5) + (c * sin(u(2))) / (m * x(3));
    J(3, 5) = x(4) + (a * c * sin(u(2))) / (m * x(3));
    J(4, 3) = -x(5)
              + (c * 1.0 / pow(x(3), 2) * (x(4) - b * x(5))
                 + c * 1.0 / pow(x(3), 2) * cos(u(2)) * (x(4) + a * x(5)))
                    / m;
    J(4, 4) = -(c / x(3) + (c * cos(u(2))) / x(3)) / m;
    J(4, 5) = -x(3) + ((b * c) / x(3) - (a * c * cos(u(2))) / x(3)) / m;
    J(5, 3) = (a * c * 1.0 / pow(x(3), 2) * (x(4) + a * x(5))
               - b * c * 1.0 / pow(x(3), 2) * (x(4) - b * x(5)))
              / I;
    J(5, 4) = -((a * c) / x(3) - (b * c) / x(3)) / I;
    J(5, 5) = -(((a * a) * c) / x(3) + ((b * b) * c) / x(3)) / I;
    return J;
}

VectorXd vehicle_observation_function(const VectorXd& x) { return x.tail(3); }

MatrixXd vehicle_observation_jacobian(const VectorXd& x)
{
    MatrixXd H = MatrixXd::Zero(3, 6);
    H(0, 3) = 1.0;
    H(1, 4) = 1.0;
    H(2, 5) = 1.0;
    return H;
}

VectorXd vehicle_cloning_state_transition_function(const VectorXd& x, const VectorXd& u)
{
    double a = 1.006;  // CG to front axle
    double b = 1.534;  // CG to rear axle

    // Construct F
    VectorXd f(6);
    f << u(0) * cos(x(2)) - u(1) * sin(x(2)), u(1) * cos(x(2)) + u(0) * sin(x(2)), u(2), 0.0, 0.0,
        0.0;

    return f;
};

MatrixXd vehicle_cloning_state_transition_jacobian(const VectorXd& x, const VectorXd& u)
{
    // Construct J
    MatrixXd J = MatrixXd::Zero(6, 6);
    J(0, 2) = -u(1) * cos(x(2)) - u(0) * sin(x(2));
    J(1, 2) = u(0) * cos(x(2)) - u(1) * sin(x(2));
    return J;
};

VectorXd vehicle_cloning_observation_function(const VectorXd& x)
{
    VectorXd h(3);
    h << cos(x(5)) * (x(0) - x(3)) + sin(x(5)) * (x(1) - x(4)),
        cos(x(5)) * (x(1) - x(4)) - sin(x(5)) * (x(0) - x(3)), x(2) - x(5);

    return h;
}

MatrixXd vehicle_cloning_observation_jacobian(const VectorXd& x)
{
    MatrixXd H = MatrixXd::Zero(3, 6);
    H(0, 0) = cos(x(5));
    H(0, 1) = sin(x(5));
    H(0, 3) = -cos(x(5));
    H(0, 4) = -sin(x(5));
    H(0, 5) = cos(x(5)) * (x(1) - x(4)) - sin(x(5)) * (x(0) - x(3));
    H(1, 0) = -sin(x(5));
    H(1, 1) = cos(x(5));
    H(1, 3) = sin(x(5));
    H(1, 4) = -cos(x(5));
    H(1, 5) = -cos(x(5)) * (x(0) - x(3)) - sin(x(5)) * (x(1) - x(4));
    H(2, 2) = 1.0;
    H(2, 5) = -1.0;
    return H;
}

runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper2;

void integrate(double dt,
               const function<VectorXd(VectorXd const&, VectorXd const&)>& f,
               VectorXd& x,
               const VectorXd& u,
               normal_random_variable& v)
{
    integrate_const(stepper2, ode(f, u, &v), x, 0.0, dt, dt / 100);
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

VectorXd vecFromYAML(const YAML::Node& node)
{
    auto v = node.as<vector<double>>();
    return Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
}

VectorXd vecFromCSV(istream& file)
{
    vector<double> result;
    string line;
    getline(file, line);

    stringstream lineStream(line);
    string cell;

    while (getline(lineStream, cell, ',')) result.push_back(stod(cell));
    return Eigen::Map<VectorXd>(result.data(), result.size());
}
