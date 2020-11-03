#include <cassert>
#include <cmath>
#include <cstdlib>  // for strtol()
#include <filesystem>
#include <fstream>
#include <iostream>
#include <utility>

#define eigen_assert(X)                         \
    do                                          \
    {                                           \
        if (!(X)) throw std::runtime_error(#X); \
    } while (false);

#include "functions/functions.h"
#include "library/library.h"
#include "third-party/multivar_noise.h"
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

    int N = (int)configCase["N"].as<double>();  /// state size
    double dt = configCase["dt"].as<double>();  /// time step
    double T = configCase["T"].as<double>();    /// total simulation time

    VectorXd x0 = vecFromYAML(configCase["x0"]);              /// Starting position
    VectorXd u = vecFromYAML(configCase["u"]);                /// Input
    MatrixXd Q = vecFromYAML(configCase["Qd"]).asDiagonal();  /// Q matrix (state transition noise)
    Q *= configCase["Qm"].as<double>();
    /// Make Q into positive semi-definite matrix
    Q = Q.transpose().eval() * Q;

    MatrixXd C =
        MatrixXd::Zero(u.size(), u.size());  /// Input noise (only used for stochastic cloning)
    if (configCase["Cd"].IsDefined())
    {
        C = vecFromYAML(configCase["Cd"]).asDiagonal();
        C *= configCase["Cm"].as<double>();
        /// Make C into positive semi-definite matrix
        C = C.transpose().eval() * C;
    }

    int nObs = configCase["obs"].size();
    vector<Observation> observations(nObs);  /// Observations
    for (int i = 0; i < nObs; i++)
    {
        MatrixXd R = vecFromYAML(configCase["obs"][i]["Rd"]).asDiagonal();  /// Observation noise
        R *= configCase["obs"][i]["Rm"].as<double>();
        /// Make R into positive semi-definite matrix
        R = R.transpose().eval() * R;
        int every_X = 1;  /// Observation occurs every 'every_X' steps
        if (configCase["obs"][i]["every_X"].IsDefined())
            every_X = (int)configCase["obs"][i]["every_X"].as<double>();
        /// Call the constructor, this will also generate the noise generator from R
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
    assert(GT_ofile.good());

    ifstream GT_ifile;
    bool GT_from_file = configCase["GT_from_file"].IsDefined();
    if (GT_from_file) GT_ifile.open(configCase["GT_from_file"].as<string>());
    assert(GT_ifile.good());

    ofstream U_ofile;
    bool U_to_file = configCase["U_to_file"].IsDefined();
    if (U_to_file) U_ofile.open(configCase["U_to_file"].as<string>(), ios::trunc);
    assert(U_ofile.good());

    ifstream U_ifile;
    bool U_from_file = configCase["U_from_file"].IsDefined();
    if (U_from_file) U_ifile.open(configCase["U_from_file"].as<string>());
    bool U_is_diff =
        U_from_file and (configCase["U_from_file"].as<string>().find("IMU") != string::npos);
    assert(U_ifile.good());

    ifstream Time_ifile;
    bool Time_from_file = configCase["Time_from_file"].IsDefined();
    if (Time_from_file) Time_ifile.open(configCase["Time_from_file"].as<string>());
    assert(Time_ifile.good());

    ofstream Kalman_ofile;
    bool Kalman_to_file = configCase["Kalman_to_file"].IsDefined();
    if (Kalman_to_file) Kalman_ofile.open(configCase["Kalman_to_file"].as<string>(), ios::trunc);
    assert(Kalman_ofile.good());

    for (unsigned int i = 0; i < observations.size(); i++)
    {
        observations[i].to_file = configCase["obs"][i]["to_file"].IsDefined();
        if (observations[i].to_file)
            observations[i].ofile.open(configCase["obs"][i]["to_file"].as<string>(), ios::trunc);
        assert(observations[i].ofile.good());
        observations[i].from_file = configCase["obs"][i]["from_file"].IsDefined();
        if (observations[i].from_file)
            observations[i].ifile.open(configCase["obs"][i]["from_file"].as<string>());
        assert(observations[i].ifile.good());
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
    bool vehicle_case = false;
    bool velocity_measurements = false;
    bool stochastic_cloning = false;
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
            vehicle_case = true;
            velocity_measurements = true;
            break;
        case 4:  /// Vehicle stochastic cloning
            ekf = ContinuousEKF(N, x0, MatrixXd::Zero(N, N));
            ekf.stochastic_cloning(0, 3, 3);
            ekf.setConstantDt(dt);
            /// In stochastic cloning the model has a different equation than the kalman filter
            /// The ground truth uses the dynamic equations
            /// The kalman filter uses kinematic equations and gets velocities from the ground truth
            /// + noise
            ekf.state_transition_function = vehicle_cloning_state_transition_function;
            f = vehicle_state_transition_function;
            ekf.state_transition_jacobian = vehicle_cloning_state_transition_jacobian;
            ekf.observation_function = vehicle_cloning_observation_function;
            ekf.observation_jacobian = vehicle_cloning_observation_jacobian;
            max_steering_angle = configCase["max_steering_angle"].as<double>();
            vehicle_case = true;
            stochastic_cloning = true;
            ekf.noise_is_input_noise = true;
            break;
        case 5:  /// Vehicle3 Smoothing
            ekf = ContinuousEKF(N, x0);
            ekf.state_transition_function = vehicle3_state_transition_function;
            f = vehicle3_state_transition_function;
            ekf.state_transition_jacobian = vehicle3_state_transition_jacobian;
            ekf.setObservationIsAngle(
                (Array<bool, Dynamic, 1>(6) << false, false, false, true, true, true).finished());
            break;
        case 6:  /// Vehicle3 Cloning
            if (GT_from_file)
            {
                x0 = vecFromCSV(GT_ifile);
                for (Observation& observation : observations)
                {
                    observation.truth_prev = x0;
                    observation.x = x0;
                }
            }

            ekf = ContinuousEKF(N, x0, MatrixXd::Zero(N, N));
            ekf.stochastic_cloning(0, N / 2, N / 2);
            ekf.setConstantDt(dt);

            ekf.state_transition_function = vehicle3_cloning_state_transition_function;
            f = vehicle3_state_transition_function;
            ekf.state_transition_jacobian = vehicle3_cloning_state_transition_jacobian;
            ekf.state_transition_input_jacobian = vehicle3_cloning_state_transition_input_jacobian;
            ekf.observation_function = vehicle3_cloning_observation_function;
            ekf.observation_jacobian = vehicle3_cloning_observation_jacobian;
            vehicle_case = true;
            stochastic_cloning = true;
            ekf.noise_is_input_noise = true;
            break;
        default: exit(1);
    }
    /// Setup integrator dt
    if (configCase["integration_dt"].IsDefined())
        ekf.setIntegrationDt(configCase["integration_dt"].as<double>());

    /// Initialize noise generators for v and w
    normal_random_variable v{Q};
    normal_random_variable zero_noise{MatrixXd::Zero(N, N)};
    normal_random_variable c{C};

    /// Start simulation
    VectorXd ground_truth = x0;
    VectorXd model_only = x0;

    double output_time = 0;
    double time = 0;
    int i = 0;
    while (time < T)
    {
        if (Time_from_file)
        {
            time = vecFromCSV(Time_ifile)(0);
        }
        else
        {
            time += dt;
        }
        i++;

        if (U_from_file)
        {
            u = vecFromCSV(U_ifile);
            if (U_is_diff)  /// If the input comes from an observation (rel pose diff) divide it by
                /// dt to get rel vel
                u /= dt;
        }
        else if (vehicle_case and i % 10 == 0 and u.size() >= 2)
        {  /// On the vehicle testcase switch steering angle randomly every 10 steps
            u(2) += fRand(-1., 1.) * 5 * EIGEN_PI / 180;  // NOLINT(cert-msc30-c,cert-msc50-cpp)
            u(2) = clamp(u(2), -max_steering_angle, max_steering_angle);
        }
        /// Simulate the movement of ground truth with added noise.
        if (GT_from_file)
            ground_truth = vecFromCSV(GT_ifile);
        else
            integrate(dt, f, ground_truth, u, v);

        VectorXd kalman_input = u;
        MatrixXd kalman_model_uncertainty = Q;
        if (vehicle_case and stochastic_cloning)
        {
            /// If we are on stochastic cloning, we only have
            kalman_model_uncertainty = C;
            if (not U_from_file)
            {
                kalman_input = ground_truth.tail(N / 2) + c();
            }
            integrate(dt, ekf.state_transition_function, model_only, kalman_input, zero_noise);
        }
        /// Call kalman prediction step.
        if (!only_ground_truth)
        {
            if (Time_from_file)
            {
                while (Kalman_to_file and time - output_time >= dt)
                {
                    output_time += dt;
                    ekf.setExternalTime(output_time);
                    ekf.predict(kalman_input, kalman_model_uncertainty);
                    Kalman_ofile << ekf.getx().format(csv) << endl;
                }
                ekf.setExternalTime(time);
            }
            ekf.predict(kalman_input, kalman_model_uncertainty);
        }

        if (isnan(ekf.getP()(1, 1))) cout << "isnan after predict" << endl;
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

                /// The z is stored as a rel pose diff in the files (to be compatible with
                /// stochastic cloning)
                if (velocity_measurements) z /= dt * observation.every_x;
            }
            else
            {
                /// Get observation from ground truth
                if (not vehicle_case)
                    z = ekf.observation_function(ground_truth);
                else if (N == 6)
                {
                    /// On the vehicle 2 test case, get absolute pose difference
                    VectorXd abs_pose_diff =
                        (ground_truth.head(N / 2) - observation.truth_prev.head(N / 2));
                    /// And then transform it to relative pose difference (z)
                    double th = observation.truth_prev(2);
                    z = VectorXd(3);
                    z << abs_pose_diff(0) * cos(th) + abs_pose_diff(1) * sin(th),  // Vn
                        -abs_pose_diff(0) * sin(th) + abs_pose_diff(1) * cos(th),  // Ve
                        abs_pose_diff(2);                                          // th_dot

                    if (velocity_measurements) z /= dt * observation.every_x;
                    observation.truth_prev = ground_truth;
                }
                else
                {
                    assert(N == 12);
                    VectorXd two_poses(12);
                    two_poses << ground_truth.head(N / 2), observation.truth_prev.head(N / 2);
                    z = ekf.observation_function(two_poses);
                    if (velocity_measurements) z /= dt * observation.every_x;
                    observation.truth_prev = ground_truth;
                }
                /// If not using stochastic cloning divide by dt to get a velocity
                /// Add noise to the observation
                z += observation.w();
            }

            /// Call kalman update step
            if (!only_ground_truth) ekf.update(z, observation.R);
            if (isnan(ekf.getP()(1, 1))) cout << "isnan after update" << endl;

            /// The z is stored as a rel pose diff in the files (to be compatible with stochastic
            /// cloning)
            if (velocity_measurements) z *= dt * observation.every_x;

            if (stochastic_cloning) ekf.stochastic_cloning(0, N / 2, N / 2);

            /// Add z to the observation
            if (not vehicle_case)
                observation.x.head(N / 2) = z;
            else if (N == 6)
            {
                /// On the vehicle test case, get absolute pose difference from rel pose diff (z)
                double ps = observation.x(2);
                VectorXd abs_pose_diff(3);
                abs_pose_diff << z(0) * cos(ps) - z(1) * sin(ps), z(0) * sin(ps) + z(1) * cos(ps),
                    z(2);
                observation.x.head(N / 2) += abs_pose_diff;
            }
            else
            {
                double ps = observation.x(3);
                double th = observation.x(4);
                double ph = observation.x(5);
                VectorXd abs_pose_diff(6);
                abs_pose_diff << -z(1) * (cos(ph) * sin(ps) - cos(ps) * sin(ph) * sin(th))
                                     + z(2) * (sin(ph) * sin(ps) + cos(ph) * cos(ps) * sin(th))
                                     + z(0) * cos(ps) * cos(th),
                    z(1) * (cos(ph) * cos(ps) + sin(ph) * sin(ps) * sin(th))
                        - z(2) * (cos(ps) * sin(ph) - cos(ph) * sin(ps) * sin(th))
                        + z(0) * cos(th) * sin(ps),
                    -z(0) * sin(th) + z(2) * cos(ph) * cos(th) + z(1) * cos(th) * sin(ph),
                    (z(5) * cos(ph)) / cos(th) + (z(4) * sin(ph)) / cos(th),
                    z(4) * cos(ph) - z(5) * sin(ph),
                    z(3) + z(5) * cos(ph) * tan(th) + z(4) * sin(ph) * tan(th);
                observation.x.head(N / 2) += abs_pose_diff;
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
        if (i % 1000 == 0) cout << i << endl;
        ofile << ground_truth.head(2).format(csv) << ",";
        ofile << nObs + ((int)(vehicle_case and stochastic_cloning)) << ",";
        if (vehicle_case and stochastic_cloning) ofile << model_only.head(2).format(csv) << ",";
        for (int iObs = 0; iObs < nObs; iObs++)
            ofile << observations[iObs].x.head(2).format(csv) << ",";
        ofile << ekf.getx().head(2).format(csv) << ",";
        ofile << ekf.getP().topLeftCorner(2, 2).format(csv) << endl;
        if (GT_to_file) GT_ofile << ground_truth.format(csv) << endl;
        if (U_to_file) U_ofile << u.format(csv) << endl;
        if (Kalman_to_file and not Time_from_file) Kalman_ofile << ekf.getx().format(csv) << endl;
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
