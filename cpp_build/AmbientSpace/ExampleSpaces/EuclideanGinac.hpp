#include <vector>
#include <ginac/ginac.h>
#include "../../includes/eigen-3.4.0/Eigen/Dense"
#include "../../Computation/State.h"
#include "../../Computation/DState.h"

using namespace GiNaC;

// Geometry Class Declaration using GiNaC

symbol a("a"), b("b"), g("g");
symbol a1("a1"), b1("b1"), g1("g1");
symbol a2("a2"), b2("b2"), g2("g2");

// Metric

matrix metric = {
    {1.,0.,0.},
    {0.,pow(sinh(a),2),0.},
    {0.,0.,pow(sinh(a),2)*pow(sin(b),2)}
};

matrix inverse_metric = metric.inverse();

ex christFunc(int i, int j, int k){
    ex christ = .5*inverse_metric[i,0];
};

// class GeometryFunctions
// {
//     public:

//         // Properties
//         std::function<Eigen::Matrix<double, 6, 6>(State)> vertex_jac;

//         // Only needed if coupling potential present
//         std::function<double(State, State)> d12;

//         std::function<double(State, State, Eigen::Vector3d)> dtd12;

//         std::function<double(State, State)> da1d12;
//         std::function<double(State, State)> db1d12;
//         std::function<double(State, State)> dg1d12;

//         std::function<double(State, State, Eigen::Vector3d, Eigen::Matrix3d)> ddtd12;

//         std::function<double(State, State)> da1da1d12;
//         std::function<double(State, State)> db1da1d12;
//         std::function<double(State, State)> dg1da1d12;
//         std::function<double(State, State)> da2da1d12;
//         std::function<double(State, State)> db2da1d12;
//         std::function<double(State, State)> dg2da1d12;

//         std::function<double(State, State)> db1db1d12;
//         std::function<double(State, State)> dg1db1d12;
//         std::function<double(State, State)> db2db1d12;
//         std::function<double(State, State)> dg2db1d12;

//         std::function<double(State, State)> dg1dg1d12;
//         std::function<double(State, State)> dg2dg1d12;

//         std::function<Eigen::Matrix3d (State)> dmetric_terms;

//         // m, f, k, l, d12,  da1d12
//         std::function<double(double,double,double,double,double,double)> coupling_derivative1;
//         // m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1
//         std::function<double(double,double,double,double,double,double,double,double,double)> coupling_derivative2;

//         // Only needed for running rigid body solvers
//         std::function<double(double,double)> rig_con;

//         std::function<double(double,double,double)> rig_con_tderivative1;
//         std::function<double(double,double,double,double)> rig_con_tderivative2;

//         // m, f, lam, d12,  da1d12
//         std::function<double(double,double,double,double,double)> rig_con_derivative1;
//         // m, f, lam, d12, da1d12, da2d12, da2f, da2d12da1
//         std::function<double(double,double,double,double,double,double,double,double)> rig_con_derivative2;

//         // dstate1, dstate2, d12, dtd12, dterms, da1dterms
//         std::function<double(DState,DState,double,double,Eigen::Vector3d, Eigen::Matrix3d)> rig_dtcon_derivative1_array;

//         std::function< Eigen::Matrix3d (Eigen::Vector3d) > metricTensor;
//         std::function< Eigen::Vector3d (State) > christoffel;
//         std::function< double (Eigen::Vector3d,Eigen::Vector3d)> distance;
//         GeometryData funcDict;

//         // Constructor
//         GeometryFunctions();

//         //Methods

// };
