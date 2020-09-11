#include <cmath>
#include <cstdlib>  // for strtol()
#include <filesystem>
#include <fstream>
#include <iostream>
#include <utility>
#include "functions/functions.h"
#include "library/library.h"
#include "multivar_noise.h"
#include "yaml-cpp/yaml.h"

using namespace std;

VectorXd vecFromYAML(const YAML::Node& node, int size)
{
    return Eigen::Map<Eigen::VectorXd>(node.as<vector<double>>().data(), size);
}

struct Observation
{
    Observation(const MatrixXd& r, int every_x, const VectorXd& x0)
        : R(r), every_x(every_x), x(x0), w(normal_random_variable(r)), truth_prev(x0)
    {
    }
    Observation() = default;
    MatrixXd R;
    int every_x = 1;
    VectorXd x;
    normal_random_variable w;
    VectorXd truth_prev;
};

int main(int argc, char* argv[])
{
    // Input variables
    int testCaseIndex = 3;
    if (argc >= 2) testCaseIndex = strtol(argv[1], nullptr, 10);

    bool only_ground_truth = false;
    if (argc >= 3) only_ground_truth = strtol(argv[2], nullptr, 10) == 1;

    string yamlfile = "config.yaml";
    if (argc >= 4) yamlfile = argv[3];

    // Input parameters
    if (!filesystem::exists(yamlfile)) yamlfile = "..\\" + yamlfile;
    YAML::Node config = YAML::LoadFile(yamlfile);
    YAML::Node configCase = config["testCase"][testCaseIndex - 1];

    // Setup output format and file
    IOFormat singleLine(StreamPrecision, DontAlignCols, ",\t", ";\t", "", "", "[", "]");
    IOFormat csv(FullPrecision, DontAlignCols, ",", ",", "", "", "", "");
    ofstream file;
    file.open("data.csv", ios::trunc);

    // Setup random
    int seed = (int)config["seed"].as<double>();
    if (seed == -1) seed = time(NULL);
    srand(seed);

    // Start reading params from config file
    int N = (int)configCase["N"].as<double>();
    double dt = configCase["dt"].as<double>();
    double T = configCase["T"].as<double>();
    Fusion f(N);
    f.setConstantDt(dt);
    VectorXd x0 = vecFromYAML(configCase["x0"], N);
    VectorXd u = vecFromYAML(configCase["u"], N);
    MatrixXd Q = vecFromYAML(configCase["Qd"], N).asDiagonal();
    Q *= configCase["Qm"].as<double>();

    int nObs = configCase["obs"].size();
    vector<Observation> observations(nObs);
    for (int i = 0; i < nObs; i++)
    {
        MatrixXd R = vecFromYAML(configCase["obs"][i]["Rd"], N / 2).asDiagonal();
        R *= configCase["obs"][i]["Rm"].as<double>();
        R = R.transpose() * R;
        int every_X = 1;
        if (configCase["obs"][i]["every_X"])
            every_X = (int)configCase["obs"][i]["every_X"].as<double>();
        observations[i] = Observation(R, every_X, x0);
    }
    int use_R2_every_x_steps = 0;
    double max_steering_angle = 0;
    switch (testCaseIndex)
    {
        case 1:  // constant_movement
            break;
        case 2:  // RTBP
            f.state_transition_function = RTBP_state_transition_function;
            f.state_transition_jacobian = RTBP_state_transition_jacobian;
            break;
        case 3:  // Vehicle dynamics
            f = Fusion(N, x0, MatrixXd::Zero(N, N));
            f.setConstantDt(dt);
            f.state_transition_function = vehicle_state_transition_function;
            f.state_transition_jacobian = vehicle_state_transition_jacobian;
            f.observation_function = vehicle_observation_function;
            f.observation_jacobian = vehicle_observation_jacobian;
            max_steering_angle = configCase["max_steering_angle"].as<double>();

            break;
        default: exit(1);
    }

    // Make Q and R into positive semi-definite matrices
    Q = Q.transpose() * Q;

    // Initialize noise generators for v and w
    normal_random_variable v{Q};
    normal_random_variable zero_noise{MatrixXd::Zero(N, N)};

    // Start simulation
    VectorXd ground_truth = x0;
    for (int i = 1; i <= (int)T / dt; i++)
    {
        // Simulate movement with simulated process noise and try to predict it
        if (testCaseIndex == 3 and i % 10 == 0)
        {
            u(2) += fRand(-1., 1.) * 5 * EIGEN_PI / 180;  // NOLINT(cert-msc30-c,cert-msc50-cpp)
            u(2) = clamp(u(2), -max_steering_angle, max_steering_angle);
        }
        integrate(dt, f.state_transition_function, ground_truth, u, v);
        if (!only_ground_truth)
        {
            if (isnan(f.getP()(1, 1))) cout << "Is nan before predict" << endl;
            f.predict(u, Q);
            if (isnan(f.getP()(1, 1))) cout << "Is nan after predict" << endl;
        }

        for (int iObs = 0; iObs < nObs; iObs++)
        {
            VectorXd z = f.observation_function(ground_truth);
            Observation obs = observations[iObs];
            if (i % obs.every_x != 0) continue;
            if (testCaseIndex == 3)
            {  // should be a boolean like use_deltas above
                VectorXd vars_dot =
                    (ground_truth.head(N / 2) - obs.truth_prev.head(N / 2)) / (dt * obs.every_x);
                double th = obs.truth_prev(2);
                z << vars_dot(0) * cos(th) + vars_dot(1) * sin(th),  // Vn
                    -vars_dot(0) * sin(th) + vars_dot(1) * cos(th),  // Ve
                    vars_dot(2);                                     // th_dot
                obs.truth_prev = ground_truth;
            }
            z += obs.w();
            if (testCaseIndex == 3)
            {  // should be a boolean like use_deltas above
                double th = obs.x(2);
                VectorXd vars_dot(3);
                vars_dot << z(0) * cos(th) - z(1) * sin(th), z(0) * sin(th) + z(1) * cos(th), z(2);
                obs.x.head(N / 2) += vars_dot * (dt * obs.every_x);
            }
            else
            {
                obs.x.head(N / 2) = z;
            }
            if (!only_ground_truth)
            {
                // Update kalman filter with simulated measurement noise
                f.update(z, obs.R);
                if (isnan(f.getP()(1, 1))) cout << "Is nan after update" << endl;
            }
            observations[iObs] = obs;
        }

        // Print information
        if (isnan(f.getP()(1, 1)))
        {
            cout << "isnan i=" << i << endl;
            return 0;
        }
        file << ground_truth.head(2).format(csv) << ",";
        file << nObs << ",";
        for (int iObs = 0; iObs < nObs; iObs++)
            file << observations[iObs].x.head(2).format(csv) << ",";
        file << f.getx().head(2).format(csv) << ",";
        file << f.getP().topLeftCorner(2, 2).format(csv) << endl;
    }
    file.close();
    return 0;
}
