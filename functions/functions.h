//
// Created by Inigo on 23/08/2020.
//

#ifndef FUSION_VISO2_IMU_FUNCTIONS_H
#define FUSION_VISO2_IMU_FUNCTIONS_H

#include <Eigen/Dense>                                            //Matrices
#include <boost/numeric/odeint.hpp>                               //Integrate
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>  //to use Matrices inside odeint
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <ctime>       //Time
#include <functional>  //Functions as parameters
#include "../third-party/multivar_noise.h"
#include "yaml-cpp/yaml.h"

using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;

VectorXd RTBP_state_transition_function(const VectorXd& x, const VectorXd& u);

MatrixXd RTBP_state_transition_jacobian(const VectorXd& x, const VectorXd& u);

VectorXd vehicle_state_transition_function(const VectorXd& x, const VectorXd& u);

MatrixXd vehicle_state_transition_jacobian(const VectorXd& x, const VectorXd& u);

VectorXd vehicle_observation_function(const VectorXd& x);

MatrixXd vehicle_observation_jacobian(const VectorXd& x);

VectorXd vehicle_cloning_state_transition_function(const VectorXd& x, const VectorXd& u);

MatrixXd vehicle_cloning_state_transition_jacobian(const VectorXd& x, const VectorXd& u);

VectorXd vehicle_cloning_observation_function(const VectorXd& x);

MatrixXd vehicle_cloning_observation_jacobian(const VectorXd& x);

VectorXd vehicle3_state_transition_function(const VectorXd& x, const VectorXd& u);

MatrixXd vehicle3_state_transition_jacobian(const VectorXd& x, const VectorXd& u);

VectorXd vehicle3_cloning_state_transition_function(const VectorXd& x, const VectorXd& u);

MatrixXd vehicle3_cloning_state_transition_jacobian(const VectorXd& x, const VectorXd& u);

MatrixXd vehicle3_cloning_state_transition_input_jacobian(const VectorXd& x);

VectorXd vehicle3_cloning_observation_function(const VectorXd& x);

MatrixXd vehicle3_cloning_observation_jacobian(const VectorXd& x);

typedef VectorXd state_type;

struct ode
{
    function<VectorXd(VectorXd const&, VectorXd const&)> f;
    VectorXd u;
    normal_random_variable* v;

    ode(function<VectorXd(VectorXd const&, VectorXd const&)> f,
        VectorXd u,
        normal_random_variable* v)
        : f(move(f)), u(move(u)), v(v)
    {
    }

    void operator()(state_type const& x, state_type& dxdt, double t) const {
        dxdt = f(x, u);
        dxdt += v->operator()();
    }
};
void integrate(double dt,
               const function<VectorXd(VectorXd const&, VectorXd const&)>& f,
               VectorXd& x,
               const VectorXd& u,
               normal_random_variable& v);
double fRand(double fMin, double fMax);

VectorXd vecFromYAML(const YAML::Node& node);

VectorXd vecFromCSV(istream& file, bool peek = false);
#endif  // FUSION_VISO2_IMU_FUNCTIONS_H
