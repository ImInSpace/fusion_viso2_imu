#ifndef FUSION_VISO2_IMU_LIBRARY_H
#define FUSION_VISO2_IMU_LIBRARY_H

#include <Eigen/Dense>  //Matrices
#include <ctime>        //Time
#include <functional>   //Functions as parameters


using namespace std;
using namespace Eigen;


class Fusion {
private:
    int N;              //State size
    clock_t t;          //Current time
    VectorXd x;         //Current state
    MatrixXd P;         //Current covariance

    //Variables to use constant dt instead of a
    bool use_constant_dt = false;
    double constant_dt = 0;

    static MatrixXd Z(int N) { return MatrixXd::Zero(N, N); } //Zero Matrix
    static MatrixXd I(int N) { return MatrixXd::Identity(N, N); } //Identity Matrix
public:
    Fusion(int N, VectorXd x, MatrixXd P) :
            N(N), t(clock()), x(move(x)), P(move(P)) {};

    explicit Fusion(int N, VectorXd x) : Fusion(N, move(x), 1000 * MatrixXd::Identity(N, N)) {};

    explicit Fusion(int N) : Fusion(N, VectorXd::Zero(N)) {};

    ~Fusion() = default;

    double predict(const VectorXd &u, const MatrixXd &Q);
    // Does the predict section of the kalman filter
    // u: input
    // Q: state transition noise
    // returns dt: time between this and last predict

    void update(const VectorXd &z, const MatrixXd &R);
    // Does the update section of the kalman filter
    // z: observation
    // R: observation noise

    [[nodiscard]] int getN() const { return N; }

    [[nodiscard]] VectorXd getx() const { return x; }

    [[nodiscard]] MatrixXd getP() const { return P; }

    void setConstantDt(float constantDt) {
        constant_dt = constantDt;
        use_constant_dt = true;
    }

    double getDt(){
        if (use_constant_dt)
            return constant_dt;
        clock_t end = clock();
        double dt = static_cast<double>(end - t) / CLOCKS_PER_SEC;
        t = end;
        return dt;
    }
};

#endif //FUSION_VISO2_IMU_LIBRARY_H