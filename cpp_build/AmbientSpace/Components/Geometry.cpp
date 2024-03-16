#include <vector>
#include <typeinfo>
#include <variant>
#include "Geometry.h"
#include "../Computation/State.h"
#include "../Computation/DState.h"

// AmbientSpace Class Declaration

struct GeometryData {
    std::function<double(State)> vertex_jac;

    std::function<double(State, State)> d12;

    std::function<double(State, State, std::vector<double>)> dtd12;

    std::function<double(State, State)> da1d12;
    std::function<double(State, State)> db1d12;
    std::function<double(State, State)> dg1d12;

    std::function<double(State, State, std::vector<double>, std::vector<std::vector<double>>)> ddtd12;

    std::function<double(State, State)> da1da1d12;
    std::function<double(State, State)> db1da1d12;
    std::function<double(State, State)> dg1da1d12;
    std::function<double(State, State)> da2da1d12;
    std::function<double(State, State)> db2da1d12;
    std::function<double(State, State)> dg2da1d12;

    std::function<double(State, State)> db1db1d12;
    std::function<double(State, State)> dg1db1d12;
    std::function<double(State, State)> db2db1d12;
    std::function<double(State, State)> dg2db1d12;

    std::function<double(State, State)> dg1dg1d12;
    std::function<double(State, State)> dg2dg1d12;

    std::function<std::vector<std::vector<double>>(State)> dmetric_terms;

    // m, f, k, l, d12,  da1d12
    std::function<double(double,double,double,double,double,double)> coupling_derivative1;
    // m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1
    std::function<double(double,double,double,double,double,double,double,double,double)> coupling_derivative2;

    std::function<double(double,double)> rig_con;

    std::function<double(double,double,double)> rig_con_tderivative1;
    std::function<double(double,double,double,double)> rig_con_tderivative2;

    // m, f, lam, d12,  da1d12
    std::function<double(double,double,double,double,double)> rig_con_derivative1;
    // m, f, lam, d12, da1d12, da2d12, da2f, da2d12da1
    std::function<double(double,double,double,double,double,double,double,double)> rig_con_derivative2;

    // dstate1, dstate2, d12, dtd12, dterms, da1dterms
    std::function<double(DState,DState,double,double,std::vector<double>, std::vector<std::vector<double>>)> rig_dtcon_derivative1_array;
};


Geometry::Geometry(
            std::function<std::vector<std::vector<double>>(std::vector<double>)> metricTensor, 
            std::function<std::vector<double>(State)> christoffel, 
            std::function<double(std::vector<double>,std::vector<double>)> distance, 
            GeometryData funcDict
            )
{
    this->metricTensor = metricTensor;
    this->christoffel = christoffel;
    this->distance = distance;
    this->funcDict = funcDict;
}

//Methods

DState Geometry::covariantAcceleration(State state)
{
    std::vector<double> acc = this->christoffel(state);
    return DState(state.vel,acc);
}

double dot(State state1, State state2);

std::vector<double> tangentBasis(std::vector<double> pos);

State gradient(std::function<double(std::vector<double>)> fn, std::vector<double> pos);



    // def covariantAcceleration(self, state):
    //     vel = state.vel
    //     acc = self.christoffel(state)
    //     return dState(vel,acc)


    // def dot(self, state1, state2):
    //     mat = self.metricTensor(state1.pos)
    //     v1 = state1.vel.copy()
    //     v2 = state2.vel.copy()
    //     # Apply this to the second vector
    //     gv2 = mat @ v2
    //     # Compute the dot product
    //     return v1.dot(gv2)


    // # Get a basis for the tangent space at a point ( in coordinates )
    // # @@@@ Consider changing how we want to do this, depending on which implementation of the gradient we want below:
    // # @@@@ right now, this is the coordinate basis, and we use the metric tensor
    // # @@@@ in the gradient: could instead do Gram-Schmidt here then calculate
    // # @@@@ gradient as differential directly in that basis.
    // def tangentBasis(self,pos):
    //     b1 = State(pos, np.array([1,0,0]))
    //     b2 = State(pos, np.array([0,1,0]))
    //     b3 = State(pos, np.array([0,0,1]))
    //     return np.array([b1,b2,b3])


    // # //WARNING: IF THE COORDINATE SYSTEM IS SINGULAR: THIS COMPUTATION IS BAD AT THAT POINT!
    // # //NEED GOOD COORDINATES.....
    // # //get the gradient of a function fn at a position pos
    // def gradient(self, fn, pos):

    //     basis = self.tangentBasis(pos)
    //     differential = State(pos, np.array([0,0,0]))

    //     #//add them all up:
    //     df0 = basis[0].differentiate(fn)
    //     b0 = basis[0].clone().multiplyScalar(df0)
    //     differential.add(b0)

    //     df1 = basis[1].differentiate(fn)
    //     b1 = basis[1].clone().multiplyScalar(df1)
    //     differential.add(b1)

    //     df2 = basis[2].differentiate(fn)
    //     b2 = basis[2].clone().multiplyScalar(df2)
    //     differential.add(b2)

    //     #//now the differential needs to be converted from a covector to a vector
    //     #//using the inverse metric:
    //     metric = self.metricTensor(pos)
    //     if(abs(np.linalg.det(metric))<0.00001):
    //         print('Warning! Metric Tensor Near Singular')
    //         print(pos)
    //         print(metric)

    //     invMetric = np.linalg.inv(metric.copy())
    //     differential.vel = invMetric @ differential.vel.copy()

    //     return differential