#ifndef FUSION_VISO2_IMU_LIBRARY_H
#define FUSION_VISO2_IMU_LIBRARY_H

#include <Eigen/Dense>  //Matrices
#include <ctime>        //Time
#include <functional>   //Functions as parameters
#include <boost/numeric/odeint.hpp> //Integrate
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp> //to use Matrices inside odeint


using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;


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

    //Functions to set and get dt, notice how if you don't set dt, you will use the clock
    void setConstantDt(double constantDt) {
        constant_dt = constantDt;
        use_constant_dt = true;
    }

    double getDt() {
        if (use_constant_dt)
            return constant_dt;
        clock_t end = clock();
        double dt = static_cast<double>(end - t) / CLOCKS_PER_SEC;
        t = end;
        return dt;
    }

    //Functions for getting state transition and observation (they need to be functions as they could change with x and u)
    function<VectorXd(VectorXd const &, VectorXd const &)> state_transition_function = [](VectorXd _x, VectorXd _u) {
        //Default f is x_dot=u+[x3;x4;0;0] (Constant velocity)
        int _N = _x.size();
        VectorXd f = _u;
        f.head(_N / 2) += _x.tail(_N / 2);
        return f;
    };

    function<MatrixXd(VectorXd const &, VectorXd const &)> state_transition_jacobian = [](VectorXd _x, VectorXd _u) {
        //Jacobian of above state transition function
        int _N = _x.size();
        MatrixXd F(_N, _N);
        F << Z(_N / 2), I(_N / 2),
                Z(_N / 2), Z(_N / 2);
        return F;
    };

    function<VectorXd(VectorXd)> observation_function = [](VectorXd _x) {
        //Only observe position, not velocity
        int _N = _x.size();
        VectorXd h = _x.head(_N/2);
        return h;
    };

    function<MatrixXd(VectorXd const &)> observation_jacobian = [](VectorXd _x) {
        //Jacobian of above observation function
        int _N = _x.size();
        MatrixXd H(_N / 2, _N);
        H << I(_N / 2), Z(_N / 2);
        return H;
    };

    //Things needed for integration
    typedef MatrixXd state_type;
    runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper;
    struct ode {
        Fusion *f;
        VectorXd u;
        MatrixXd Q;

        ode(Fusion *f, VectorXd u, MatrixXd Q) : f(f), u(move(u)), Q(move(Q)) {}

        void operator()(state_type const &pair, state_type &dpairdt, double t) const {
            int N = pair.rows();
            dpairdt = state_type(N, N + 1);
            dpairdt.col(0) = f->state_transition_function(pair.col(0), u);
            MatrixXd F = f->state_transition_jacobian(pair.col(0), u);
            dpairdt.rightCols(N) = F * pair.rightCols(N) + pair.rightCols(N) * F.transpose() + Q;
        }
    };

};

#endif //FUSION_VISO2_IMU_LIBRARY_H