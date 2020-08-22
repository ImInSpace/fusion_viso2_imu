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
protected:

    /** State size **/
    int N;

    /** Current time **/
    clock_t t;

    /** Current state **/
    VectorXd x;

    /** Current covariance **/
    MatrixXd P;

    /** Variables to use constant dt instead of getting it from clock() **/
    bool use_constant_dt = false;
    double constant_dt = 0;

    /** Returns time to be predicted in a predict step. **/
    double getDt() {
        if (use_constant_dt)
            return constant_dt;
        clock_t end = clock();
        double dt = static_cast<double>(end - t) / CLOCKS_PER_SEC;
        t = end;
        return dt;
    }


    /** Helper functions **/
    static MatrixXd Z(int N) { return MatrixXd::Zero(N, N); }
    static MatrixXd I(int N) { return MatrixXd::Identity(N, N); }
public:

    /** Constructor for Fusion
     * \param N size of the state vector.
     * \param x initial state of the kalman prediction. Default is zero.
     * \param P initial covariance matrix for the kalman prediction. Default is 1000*Identity.
     */
    Fusion(int N, VectorXd x, MatrixXd P) :
            N(N), t(clock()), x(move(x)), P(move(P)) {};

    explicit Fusion(int N, VectorXd x) : Fusion(N, move(x), 1000 * MatrixXd::Identity(N, N)) {};

    explicit Fusion(int N) : Fusion(N, VectorXd::Zero(N)) {};

    ~Fusion() = default;

    /** Does the predict section of the kalman filter.
     * \param u input.
     * \param Q state transition noise.
     * \return dt distance of time predicted
     */
    double predict(const VectorXd &u, const MatrixXd &Q);

    /** Does the update section of the kalman filter
     * @param z observation
     * @param R observation noise
     */
    void update(const VectorXd &z, const MatrixXd &R);

    [[nodiscard]] int getN() const { return N; }

    [[nodiscard]] VectorXd getx() const { return x; }

    [[nodiscard]] MatrixXd getP() const { return P; }

    /** Set dt (the time predicted by each predict step to a constant value. If not used, it will use the clock
     * @param constantDt value of dt to be used.
     */
    void setConstantDt(double constantDt) {
        constant_dt = constantDt;
        use_constant_dt = true;
    }

    /** State transition function f: x'=f(x,u).
     * Default f is x_dot=u+[x3;x4;0;0] (Constant velocity)
     */
    function<VectorXd(VectorXd const &, VectorXd const &)> state_transition_function = [](VectorXd _x, VectorXd _u) {
        int _N = _x.size();
        VectorXd f = _u;
        f.head(_N / 2) += _x.tail(_N / 2);
        return f;
    };

    /** Jacobian of state_transition_function **/
    function<MatrixXd(VectorXd const &, VectorXd const &)> state_transition_jacobian = [](VectorXd _x, VectorXd _u) {
        int _N = _x.size();
        MatrixXd F(_N, _N);
        F << Z(_N / 2), I(_N / 2),
                Z(_N / 2), Z(_N / 2);
        return F;
    };

    /** observation function h: z=h(x).
     * Default is to only observe first half of x
     */
    function<VectorXd(VectorXd)> observation_function = [](VectorXd _x) {
        int _N = _x.size();
        VectorXd h = _x.head(_N/2);
        return h;
    };

    /** Jacobian of the observation function **/
    function<MatrixXd(VectorXd const &)> observation_jacobian = [](VectorXd _x) {
        //Jacobian of above observation function
        int _N = _x.size();
        MatrixXd H(_N / 2, _N);
        H << I(_N / 2), Z(_N / 2);
        return H;
    };

    /** Things needed for integration **/
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