#ifndef FUSION_VISO2_IMU_LIBRARY_H
#define FUSION_VISO2_IMU_LIBRARY_H

#include <Eigen/Dense>  //Matrices
#include <boost/config/warning_disable.hpp>
#include <boost/numeric/odeint.hpp>                               //Integrate
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>  //to use Matrices inside odeint
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <ctime>       //Time
#include <functional>  //Functions as parameters

using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;

class ContinuousEKF
{
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

    /** Variables to use external time **/
    bool use_external_t = false;
    double current_external_t = 0;  /// Where the kalman filter has integrated
    double external_t = 0;          /// The external time

    /** Number of integration steps for each predict **/
    double integration_dt = 0.001;

    /** Observation is angle **/
    typedef Array<bool, Dynamic, 1> ArrayXb;
    ArrayXb observation_is_angle;

  public:
    /** Noise is input noise instead of state transition **/
    bool noise_is_input_noise=false;

  protected:
    /** Returns time to be predicted in a predict step. **/
    double getDt()
    {

        if (use_constant_dt) return constant_dt;

        if (use_external_t)
        {
            double dt = external_t - current_external_t;
            current_external_t = external_t;
            return dt;
        }
        clock_t end = clock();
        double dt = static_cast<double>(end - t) / CLOCKS_PER_SEC;
        t = end;
        return dt;
    }

    /** Helper functions **/
    static MatrixXd Z(int N) { return MatrixXd::Zero(N, N); }

    static MatrixXd I(int N) { return MatrixXd::Identity(N, N); }

  public:
    /** Constructor for ContinuousEKF
     * \param N size of the state vector.
     * \param x initial state of the kalman prediction. Default is zero.
     * \param P initial covariance matrix for the kalman prediction. Default is 1000*Identity.
     */
    ContinuousEKF(int N, VectorXd x, MatrixXd P) : N(N), t(clock()), x(move(x)), P(move(P)){};

    explicit ContinuousEKF(int N, VectorXd x)
        : ContinuousEKF(N, move(x), 1000 * MatrixXd::Identity(N, N)){};

    explicit ContinuousEKF(int N) : ContinuousEKF(N, VectorXd::Zero(N)){};

    ~ContinuousEKF() = default;

    /** Does the predict section of the kalman filter.
     * \param u input.
     * \param Q state transition noise.
     * \return dt distance of time predicted
     */
    double predict(const VectorXd& u, const MatrixXd& Q);

    /** Does the update section of the kalman filter
     * @param z observation
     * @param R observation noise
     */
    void update(const VectorXd& z, const MatrixXd& R);

    [[nodiscard]] int getN() const { return N; }

    [[nodiscard]] VectorXd getx() const { return x; }

    [[nodiscard]] MatrixXd getP() const { return P; }

    /** Set dt (the time predicted by each predict) step to a constant value. If not used, it will
     * use the clock
     * @param constantDt value of dt to be used.
     */
    void setConstantDt(double constantDt)
    {
        constant_dt = constantDt;
        use_constant_dt = true;
    }

    /** Set externalTime (the current time) the predict step will integrate from the previous time
     * to the new time
     * @param externalTime current_time
     */
    void setExternalTime(double externalTime)
    {
        external_t = externalTime;
        use_external_t = true;
    }

    void setIntegrationDt(double integrationDt) { integration_dt = integrationDt; }

    void setObservationIsAngle(const ArrayXb& observationIsAngle)
    {
        observation_is_angle = observationIsAngle;
    }

    /** Clone a part of the kalman state, also affects the covariance
     * Basically, x(to+i)=x(from+i) for 0<=i<size (changing the covariance accordingly)
     * @param from index where to start copying from
     * @param to index where to start copying to
     * @param size number of elements to copy
     */
    void stochastic_cloning(int from, int to, int size)
    {
        x.segment(to, size) = x.segment(from, size);
        P.block(to, 0, size, N) = P.block(from, 0, size, N);
        P.block(0, to, N, size) = P.block(0, from, N, size);
    }

    typedef function<VectorXd(const VectorXd&, const VectorXd&)> State_transition_function;
    /** State transition function f: x'=f(x,u).
     * Default f is x_dot=u+[x3;x4;0;0] (Constant velocity)
     */
    State_transition_function state_transition_function = [](const VectorXd& _x,
                                                             const VectorXd& _u) {
        int _N = _x.size();
        VectorXd f = _u;
        f.head(_N / 2) += _x.tail(_N / 2);
        return f;
    };

    typedef function<MatrixXd(VectorXd const&, VectorXd const&)> State_transition_jacobian;
    /** Jacobian of state_transition_function wrt the state**/
    State_transition_jacobian state_transition_jacobian = [](const VectorXd& _x,
                                                             const VectorXd& _u) {
        int _N = _x.size();
        MatrixXd F(_N, _N);
        F << Z(_N / 2), I(_N / 2), Z(_N / 2), Z(_N / 2);
        return F;
    };

    typedef function<MatrixXd(VectorXd const&)> State_transition_input_jacobian;
    /** Jacobian of state_transition_function wrt the input**/
    State_transition_input_jacobian state_transition_input_jacobian = [](const VectorXd& _x) {
        int _N = _x.size();
        MatrixXd F(_N, _N / 2);
        F << I(_N / 2), Z(_N / 2);
        return F;
    };

    typedef function<VectorXd(VectorXd)> Observation_function;
    /** observation function h: z=h(x).
     * Default is to only observe first half of x
     */
    Observation_function observation_function = [](const VectorXd& _x) {
        int _N = _x.size();
        VectorXd h = _x.head(_N / 2);
        return h;
    };

    typedef function<MatrixXd(VectorXd const&)> Observation_jacobian;
    /** Jacobian of the observation function **/
    Observation_jacobian observation_jacobian = [](const VectorXd& _x) {
        // Jacobian of above observation function
        int _N = _x.size();
        MatrixXd H(_N / 2, _N);
        H << I(_N / 2), Z(_N / 2);
        return H;
    };

    /** Things needed for integration **/
    typedef MatrixXd state_type;
    runge_kutta_dopri5<state_type, double, state_type, double, vector_space_algebra> stepper;

    struct ode
    {
        ContinuousEKF* f;
        VectorXd u;
        MatrixXd Q;

        ode(ContinuousEKF* f, VectorXd u, MatrixXd Q) : f(f), u(move(u)), Q(move(Q)) {}

        void operator()(state_type const& pair, state_type& dpairdt, double t) const;
    };
};

#endif  // FUSION_VISO2_IMU_LIBRARY_H