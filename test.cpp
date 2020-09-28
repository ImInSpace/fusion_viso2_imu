#include <cassert>
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
    ofstream ofile;
    ifstream ifile;
    bool from_file = false;
    bool to_file = false;
};

int main(int argc, char* argv[])
{
    /// Input parameters (from argv)
    int testCaseIndex = 4;
    if (argc >= 2) testCaseIndex = strtol(argv[1], nullptr, 10);

    bool only_ground_truth = false;
    if (argc >= 3) only_ground_truth = strtol(argv[2], nullptr, 10) == 1;

    string yamlfile = "config.yaml";
    if (argc >= 4) yamlfile = argv[3];

    /// Input parameters (from YAML)
    if (!filesystem::exists(yamlfile)) filesystem::current_path("..");
    assert(filesystem::exists(yamlfile));
    YAML::Node config = YAML::LoadFile(yamlfile);
    YAML::Node configCase = config["testCase"][testCaseIndex - 1];

    int N = (int)configCase["N"].as<double>();
    double dt = configCase["dt"].as<double>();
    double T = configCase["T"].as<double>();

    VectorXd x0 = vecFromYAML(configCase["x0"]);
    VectorXd u = vecFromYAML(configCase["u"]);
    MatrixXd Q = vecFromYAML(configCase["Qd"]).asDiagonal();
    Q *= configCase["Qm"].as<double>();

    int nObs = configCase["obs"].size();
    vector<Observation> observations(nObs);
    for (int i = 0; i < nObs; i++)
    {
        MatrixXd R = vecFromYAML(configCase["obs"][i]["Rd"]).asDiagonal();
        R *= configCase["obs"][i]["Rm"].as<double>();
        R = R.transpose().eval() * R;  // Make R positive semi-definite
        int every_X = 1;
        if (configCase["obs"][i]["every_X"].IsDefined())
            every_X = (int)configCase["obs"][i]["every_X"].as<double>();
        observations[i] = Observation(R, every_X, x0);
    }

    /// Set up Eigen output formats
    IOFormat singleLine(StreamPrecision, DontAlignCols, ",\t", ";\t", "", "", "[", "]");
    IOFormat csv(FullPrecision, DontAlignCols, ",", ",", "", "", "", "");

    /// Set up input/output files
    ofstream ofile;
    ofile.open("data.csv", ios::trunc);
    ofstream GT_ofile;
    bool GT_to_file = configCase["GT_to_file"].IsDefined();
    if (GT_to_file) GT_ofile.open(configCase["GT_to_file"].as<string>(), ios::trunc);
    if (GT_to_file) assert(filesystem::exists(configCase["GT_to_file"].as<string>()));
    ifstream GT_ifile;
    bool GT_from_file = configCase["GT_from_file"].IsDefined();
    if (GT_from_file) GT_ifile.open(configCase["GT_from_file"].as<string>());
    if (GT_from_file) assert(filesystem::exists(configCase["GT_from_file"].as<string>()));
    ofstream U_ofile;
    bool U_to_file = configCase["U_to_file"].IsDefined();
    if (U_to_file) U_ofile.open(configCase["U_to_file"].as<string>(), ios::trunc);
    ifstream U_ifile;
    bool U_from_file = configCase["U_from_file"].IsDefined();
    if (U_from_file) U_ifile.open(configCase["U_from_file"].as<string>());
    bool U_is_diff =
        U_from_file and (configCase["U_from_file"].as<string>().find("IMU") != string::npos);
    for (int i = 0; i < observations.size(); i++)
    {
        observations[i].to_file = configCase["obs"][i]["to_file"].IsDefined();
        if (observations[i].to_file)
            observations[i].ofile.open(configCase["obs"][i]["to_file"].as<string>(), ios::trunc);
        observations[i].from_file = configCase["obs"][i]["from_file"].IsDefined();
        if (observations[i].from_file)
            observations[i].ifile.open(configCase["obs"][i]["from_file"].as<string>());
    }

    /// Set up random
    int seed = time(NULL);
    if (config["seed"].IsDefined()) seed = (int)config["seed"].as<double>();
    srand(seed);

    /// Add state and observation functions to the kalman library
    ContinuousEKF ekf(N);
    ekf.setConstantDt(dt);
    double max_steering_angle = 0;
    ContinuousEKF::State_transition_function f;
    switch (testCaseIndex)
    {
        case 1:  /// constant_movement
            f = ekf.state_transition_function;
            break;
        case 2:  /// RTBP
            ekf = ContinuousEKF(N, x0, MatrixXd::Zero(N, N));
            ekf.state_transition_function = RTBP_state_transition_function;
            f = RTBP_state_transition_function;
            ekf.state_transition_jacobian = RTBP_state_transition_jacobian;
            break;
        case 3:  /// Vehicle dynamics
            ekf = ContinuousEKF(N, x0, MatrixXd::Zero(N, N));
            ekf.setConstantDt(dt);
            ekf.state_transition_function = vehicle_state_transition_function;
            f = vehicle_state_transition_function;
            ekf.state_transition_jacobian = vehicle_state_transition_jacobian;
            ekf.observation_function = vehicle_observation_function;
            ekf.observation_jacobian = vehicle_observation_jacobian;
            max_steering_angle = configCase["max_steering_angle"].as<double>();
            break;
        case 4:  /// Vehicle stochastic cloning
            ekf = ContinuousEKF(N, x0, MatrixXd::Zero(N, N));
            ekf.stochastic_cloning(0, 3, 3);
            ekf.setConstantDt(dt);
            ekf.state_transition_function = vehicle_cloning_state_transition_function;
            f = vehicle_state_transition_function;
            ekf.state_transition_jacobian = vehicle_cloning_state_transition_jacobian;
            ekf.observation_function = vehicle_cloning_observation_function;
            ekf.observation_jacobian = vehicle_cloning_observation_jacobian;
            max_steering_angle = configCase["max_steering_angle"].as<double>();
            break;
        default: exit(1);
    }
    /// Setup integrator steps
    if (configCase["integration_steps"].IsDefined())
        ekf.setSteps((int)configCase["integration_steps"].as<double>());

    /// Make Q into positive semi-definite matrix
    Q = Q.transpose().eval() * Q;

    /// Initialize noise generators for v and w
    normal_random_variable v{Q};
    normal_random_variable zero_noise{MatrixXd::Zero(N, N)};

    /// Start simulation
    VectorXd ground_truth = x0;
    for (int i = 1; i <= (int)T / dt; i++)
    {
        /// On the vehicle testcase switch steering angle randomly every 10 steps
        if (U_from_file)
        {
            u = vecFromCSV(U_ifile);
            if (U_is_diff)
                u /= dt;
        }
        else if (testCaseIndex == 3 and i % 10 == 0)
        {
            u(2) += fRand(-1., 1.) * 5 * EIGEN_PI / 180;  // NOLINT(cert-msc30-c,cert-msc50-cpp)
            u(2) = clamp(u(2), -max_steering_angle, max_steering_angle);
        }
        /// Simulate the movement of ground truth with added noise.
        if (GT_from_file)
            ground_truth = vecFromCSV(GT_ifile);
        else
            integrate(dt, f, ground_truth, u, v);

        /// Call kalman prediction step.
        if (!only_ground_truth) ekf.predict(u, Q);

        /// Iterate through all observations
        for (auto& observation : observations)
        {
            /// Skip observation if it is not its turn
            if (i % observation.every_x != 0) continue;

            VectorXd z;
            /// Get observation
            if (observation.from_file)
            {
                z = vecFromCSV(observation.ifile);
                if (testCaseIndex == 3) z /= dt * observation.every_x;
            }
            else
            {
                /// Get observation from ground truth
                if (testCaseIndex != 3 and testCaseIndex != 4)
                    z = ekf.observation_function(ground_truth);
                else
                {
                    /// On the vehicle test case, get pose diff
                    VectorXd vars_dot =
                        (ground_truth.head(N / 2) - observation.truth_prev.head(N / 2));
                    /// And then transform it to local reference
                    double th = observation.truth_prev(2);
                    z = VectorXd(3);
                    z << vars_dot(0) * cos(th) + vars_dot(1) * sin(th),  // Vn
                        -vars_dot(0) * sin(th) + vars_dot(1) * cos(th),  // Ve
                        vars_dot(2);                                     // th_dot

                    if (testCaseIndex == 3) z /= dt * observation.every_x;
                    observation.truth_prev = ground_truth;
                }
                /// If not using stochastic cloning divide by dt to get a velocity
                /// Add noise to the observation
                z += observation.w();
            }

            /// Call kalman update step
            if (!only_ground_truth) ekf.update(z, observation.R);

            if (testCaseIndex == 3) z *= dt * observation.every_x;

            if (testCaseIndex == 4) ekf.stochastic_cloning(0, 3, 3);

            /// Add z to the observation
            if (testCaseIndex != 3 and testCaseIndex != 4)
                observation.x.head(N / 2) = z;
            else
            {
                /// On the vehicle test case, undo the previous transformation
                double th = observation.x(2);
                VectorXd vars_dot(3);
                vars_dot << z(0) * cos(th) - z(1) * sin(th), z(0) * sin(th) + z(1) * cos(th), z(2);
                observation.x.head(N / 2) += vars_dot;
            }

            /// Print observation to file
            if (observation.to_file) observation.ofile << z.format(csv) << endl;
        }

        /// Print information
        if (isnan(ekf.getP()(1, 1)))
        {
            cout << "isnan i=" << i << endl;
            return 0;
        }
        ofile << ground_truth.head(2).format(csv) << ",";
        ofile << nObs << ",";
        for (int iObs = 0; iObs < nObs; iObs++)
            ofile << observations[iObs].x.head(2).format(csv) << ",";
        ofile << ekf.getx().head(2).format(csv) << ",";
        ofile << ekf.getP().topLeftCorner(2, 2).format(csv) << endl;
        if (GT_to_file) GT_ofile << ground_truth.format(csv) << endl;
        if (U_to_file) U_ofile << u.format(csv) << endl;
    }
    ofile.close();
    if (GT_to_file) GT_ofile.close();
    if (GT_from_file) GT_ifile.close();
    if (U_to_file) U_ofile.close();
    if (U_from_file) U_ifile.close();
    for (auto& observation : observations)
    {
        if (observation.to_file) observation.ofile.close();
        if (observation.from_file) observation.ifile.close();
    }
    return 0;
}
