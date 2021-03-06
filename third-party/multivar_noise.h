// from
// https://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
#ifndef FUSION_VISO2_IMU_MULTIVAR_NOISE_H
#define FUSION_VISO2_IMU_MULTIVAR_NOISE_H

#include <Eigen/Dense>
#include <random>
#include <utility>

struct normal_random_variable
{
    normal_random_variable() : normal_random_variable(Eigen::MatrixXd::Identity(2, 2)) {}
    explicit normal_random_variable(Eigen::MatrixXd const& covar)
        : normal_random_variable(Eigen::VectorXd::Zero(covar.rows()), covar)
    {
    }

    normal_random_variable(Eigen::VectorXd mean, Eigen::MatrixXd const& covar)
        : mean(std::move(mean))
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
        gen = std::mt19937{static_cast<unsigned>(rand())};
    }

    Eigen::VectorXd mean;
    Eigen::MatrixXd transform;

    std::mt19937 gen;

    Eigen::VectorXd operator()()
    {
        static std::normal_distribution<> dist;
        return mean + transform * Eigen::VectorXd{mean.size()}.unaryExpr([&](auto x) {
            return dist(gen);
        });
    }
};

#endif  // FUSION_VISO2_IMU_MULTIVAR_NOISE_H