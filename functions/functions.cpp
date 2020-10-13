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

    // from matlab_functions/kinematic_functions/Vehicle.m
    VectorXd f = VectorXd::Zero(6);
    double t2 = cos(x(2));
    double t3 = sin(x(2));
    double t4 = a * x(5);
    double t5 = b * x(5);
    double t6 = 1.0 / m;
    double t7 = 1.0 / x(3);
    double t8 = t4 + x(4);
    double t9 = -t5;
    double t10 = t9 + x(4);
    double t11 = t7 * t8;
    double t12 = -t11;
    f(0) = t2 * x(3) - t3 * x(4);
    f(1) = t2 * x(4) + t3 * x(3);
    f(2) = x(5);
    f(3) = x(4) * x(5) + t6 * (u(1) + c * sin(u(2)) * (t11 - u(2)));
    f(4) = -x(3) * x(5) - t6 * (c * cos(u(2)) * (t11 - u(2)) - c * t7 * (t5 - x(4)));
    f(5) = -(a * c * (t11 - u(2)) + b * c * t7 * (t5 - x(4))) / I;

    return f;
}

MatrixXd vehicle_state_transition_jacobian(const VectorXd& x, const VectorXd& u)
{
    double m = 1292.2;     // Vehicle mass
    double I = 2380.7;     // Vehicle inertia
    double a = 1.006 / 2;  // CG to front axle
    double b = 1.534 / 2;  // CG to rear axle
    double c = 20000;      // cornering stiffness

    // from matlab_functions/kinematic_functions/Vehicle.m
    MatrixXd J = MatrixXd::Zero(6, 6);
    double t2 = cos(u(2));
    double t3 = cos(x(2));
    double t4 = sin(u(2));
    double t5 = sin(x(2));
    double t6 = a * x(5);
    double t7 = a * a;
    double t8 = b * b;
    double t9 = 1.0 / I;
    double t10 = 1.0 / m;
    double t11 = 1.0 / x(3);
    double t12 = t11 * t11;
    double t13 = t6 + x(4);
    J(0, 2) = -t3 * x(4) - t5 * x(3);
    J(0, 3) = t3;
    J(0, 4) = -t5;
    J(1, 2) = t3 * x(3) - t5 * x(4);
    J(1, 3) = t5;
    J(1, 4) = t3;
    J(2, 5) = 1.0;
    J(3, 3) = -c * t4 * t10 * t12 * t13;
    J(3, 4) = x(5) + c * t4 * t10 * t11;
    J(3, 5) = x(4) + a * c * t4 * t10 * t11;
    J(4, 3) = -x(5) + t10 * (c * t12 * (x(4) - b * x(5)) + c * t2 * t12 * t13);
    J(4, 4) = -c * t10 * t11 * (t2 + 1.0);
    J(4, 5) = -x(3) + t10 * t11 * (b * c - a * c * t2);
    J(5, 3) = c * t9 * t12 * (a * t6 + a * x(4) - b * x(4) + t8 * x(5));
    J(5, 4) = -c * t9 * t11 * (a - b);
    J(5, 5) = -c * t9 * t11 * (t7 + t8);
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
    // from matlab_functions/kinematic_functions/Vehicle_Kin_Cloning.m
    double var1 = sin(x(2));
    double var2 = cos(x(2));
    VectorXd f = VectorXd::Zero(6);
    f(0) = -var1 * u(1) + var2 * u(0);
    f(1) = var1 * u(0) + var2 * u(1);
    f(2) = u(2);

    return f;
};

MatrixXd vehicle_cloning_state_transition_jacobian(const VectorXd& x, const VectorXd& u)
{
    // from matlab_functions/kinematic_functions/Vehicle_Kin_Cloning.m
    MatrixXd J = MatrixXd::Zero(6, 6);
    double t2 = cos(x(2));
    double t3 = sin(x(2));
    J(0, 2) = -t2 * u(1) - t3 * u(0);
    J(1, 2) = t2 * u(0) - t3 * u(1);
    return J;
};

VectorXd vehicle_cloning_observation_function(const VectorXd& x)
{
    // from matlab_functions/kinematic_functions/Vehicle_Kin_Cloning.m
    VectorXd h = VectorXd::Zero(3);
    double t2 = cos(x(5));
    double t3 = sin(x(5));
    double t4 = -x(3);
    double t5 = -x(4);
    double t6 = t4 + x(0);
    double t7 = t5 + x(1);
    h(0) = t2 * t6 + t3 * t7;
    h(1) = t2 * t7 - t3 * t6;
    h(2) = x(2) - x(5);

    return h;
}

MatrixXd vehicle_cloning_observation_jacobian(const VectorXd& x)
{
    // from matlab_functions/kinematic_functions/Vehicle_Kin_Cloning.m
    MatrixXd H = MatrixXd::Zero(3, 6);
    double t2 = cos(x(5));
    double t3 = sin(x(5));
    double t4 = -x(3);
    double t5 = -x(4);
    double t6 = -t2;
    double t7 = -t3;
    double t8 = t4 + x(0);
    double t9 = t5 + x(1);
    H(0, 0) = t2;
    H(0, 1) = t3;
    H(0, 3) = t6;
    H(0, 4) = t7;
    H(0, 5) = t2 * t9 + t7 * t8;
    H(1, 0) = t7;
    H(1, 1) = t2;
    H(1, 3) = t3;
    H(1, 4) = t6;
    H(1, 5) = t6 * t8 + t7 * t9;
    H(2, 2) = 1.0;
    H(2, 5) = -1.0;
    return H;
}

VectorXd vehicle3_state_transition_function(const VectorXd& x, const VectorXd& u)
{

    // from matlab_functions/kinematic_functions/Vehicle3_Kin.m
    VectorXd f = VectorXd::Zero(12);
    double t2 = cos(x(3));
    double t3 = cos(x(4));
    double t4 = cos(x(5));
    double t5 = sin(x(3));
    double t6 = sin(x(4));
    double t7 = sin(x(5));
    double t8 = t4 * x(11);
    double t9 = t7 * x(10);
    double t10 = 1.0 / t3;
    double t11 = t8 + t9;
    f(0) = -x(7) * (t4 * t5 - t2 * t6 * t7) + x(8) * (t5 * t7 + t2 * t4 * t6) + t2 * t3 * x(6);
    f(1) = x(7) * (t2 * t4 + t5 * t6 * t7) - x(8) * (t2 * t7 - t4 * t5 * t6) + t3 * t5 * x(6);
    f(2) = -t6 * x(6) + t3 * t4 * x(8) + t3 * t7 * x(7);
    f(3) = t10 * t11;
    f(4) = t4 * x(10) - t7 * x(11);
    f(5) = x(9) + t6 * t10 * t11;

    return f;
}

MatrixXd vehicle3_state_transition_jacobian(const VectorXd& x, const VectorXd& u)
{
    // from matlab_functions/kinematic_functions/Vehicle3_Kin.m
    MatrixXd J = MatrixXd::Zero(12, 12);
    double t2 = cos(x(3));
    double t3 = cos(x(4));
    double t4 = cos(x(5));
    double t5 = sin(x(3));
    double t6 = sin(x(4));
    double t7 = sin(x(5));
    double t8 = tan(x(4));
    double t9 = t4 * x(10);
    double t10 = t4 * x(11);
    double t11 = t6 * x(6);
    double t12 = t7 * x(10);
    double t13 = t7 * x(11);
    double t14 = t2 * t4;
    double t15 = t2 * t7;
    double t16 = t4 * t5;
    double t17 = t5 * t7;
    double t18 = 1.0 / t3;
    double t19 = t3 * t4 * x(8);
    double t20 = t3 * t7 * x(7);
    double t21 = -t11;
    double t22 = -t12;
    double t23 = -t13;
    double t24 = t6 * t14;
    double t25 = t6 * t15;
    double t26 = t6 * t16;
    double t27 = t6 * t17;
    double t28 = -t25;
    double t29 = -t26;
    double t30 = t9 + t23;
    double t31 = t14 + t27;
    double t32 = t17 + t24;
    double t35 = t19 + t20 + t21;
    double t33 = t15 + t29;
    double t34 = t16 + t28;
    J(0, 3) = -t31 * x(7) + t33 * x(8) - t3 * t5 * x(6);
    J(0, 4) = t2 * t35;
    J(0, 5) = t32 * x(7) + t34 * x(8);
    J(0, 6) = t2 * t3;
    J(0, 7) = -t16 + t25;
    J(0, 8) = t32;
    J(1, 3) = t32 * x(8) - t34 * x(7) + t2 * t3 * x(6);
    J(1, 4) = t5 * t35;
    J(1, 5) = -t31 * x(8) - t33 * x(7);
    J(1, 6) = t3 * t5;
    J(1, 7) = t31;
    J(1, 8) = -t15 + t26;
    J(2, 4) = -t3 * x(6) - t4 * t6 * x(8) - t6 * t7 * x(7);
    J(2, 5) = t3 * (t4 * x(7) - t7 * x(8));
    J(2, 6) = -t6;
    J(2, 7) = t3 * t7;
    J(2, 8) = t3 * t4;
    J(3, 4) = -(t6 * (t12 - x(11) * (pow(sin(x(5) / 2.0), 2.0) * 2.0 - 1.0))) / (t6 * t6 - 1.0);
    J(3, 5) = t18 * t30;
    J(3, 10) = t7 * t18;
    J(3, 11) = t4 * t18;
    J(4, 5) = -t10 + t22;
    J(4, 10) = t4;
    J(4, 11) = -t7;
    J(5, 4) = (t18 * t18) * (t10 + t12);
    J(5, 5) = t6 * t18 * t30;
    J(5, 9) = 1.0;
    J(5, 10) = t7 * t8;
    J(5, 11) = t4 * t8;
    return J;
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
