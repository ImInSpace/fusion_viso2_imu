// from https://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
#ifndef FUSION_VISO2_IMU_MULTIVAR_NOISE_H
#define FUSION_VISO2_IMU_MULTIVAR_NOISE_H

#include <Eigen\Dense>
#include <utility>
#include <random>

struct normal_random_variable
{
    explicit normal_random_variable(Eigen::MatrixXd const& covar)
            : normal_random_variable(Eigen::VectorXd::Zero(covar.rows()), covar)
    {}

    normal_random_variable(Eigen::VectorXd  mean, Eigen::MatrixXd const& covar)
            : mean(std::move(mean))
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::VectorXd mean;
    Eigen::MatrixXd transform;

    Eigen::VectorXd operator()() const
    {
        static std::mt19937 gen{ std::random_device{}() };
        static std::normal_distribution<> dist;

        return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](auto x) { return dist(gen); });
    }
};

#endif //FUSION_VISO2_IMU_MULTIVAR_NOISE_H