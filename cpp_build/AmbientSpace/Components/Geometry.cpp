#include <vector>
#include <iostream>
#include <typeinfo>
#include <variant>
#include "../../includes/eigen-3.4.0/Eigen/Dense"
#include "Geometry.h"
#include "../../Computation/State.h"
#include "../../Computation/DState.h"

// AmbientSpace Class Declaration

Geometry::Geometry()
{

}

Geometry::Geometry(
            std::function< Eigen::Matrix3d (Eigen::Vector3d)> metricTensor, 
            std::function< Eigen::Vector3d (State)> christoffel, 
            std::function<double(Eigen::Vector3d,Eigen::Vector3d)> distance
            )
{
    this->metricTensor = metricTensor;
    this->christoffel = christoffel;
    this->distance = distance;
}

Geometry::Geometry(
            std::function< Eigen::Matrix3d (Eigen::Vector3d)> metricTensor, 
            std::function< Eigen::Vector3d (State)> christoffel, 
            std::function<double(Eigen::Vector3d,Eigen::Vector3d)> distance, 
            derivative_funcs funcDict
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
    Eigen::Vector3d acc = this->christoffel(state);
    return DState(state.vel,acc);
}

double Geometry::dot(State state1, State state2)
{
    Eigen::Matrix3d mat = this->metricTensor(state1.pos);
    
    return state1.vel.transpose() * mat * state2.vel;
}

std::vector<State> Geometry::tangentBasis(Eigen::Vector3d pos)
{
    State b1 = State(pos, Eigen::Vector3d(1.,0.,0.));
    State b2 = State(pos, Eigen::Vector3d(0.,1.,0.));
    State b3 = State(pos, Eigen::Vector3d(0.,0.,1.));
    std::vector<State> basis = {b1,b2,b3};
    return basis;
}

State Geometry::gradient(std::function<double(Eigen::Vector3d)> fn, Eigen::Vector3d pos)
{
    std::vector<State> basis = tangentBasis(pos);
    State differential = State(pos, Eigen::Vector3d(0.,0.,0.));

    // #//add them all up:
    double df0 = basis[0].differentiate(fn);
    State b0 = basis[0].clone();
    b0.multiplyScalar(df0);
    differential.add(b0);

    double df1 = basis[1].differentiate(fn);
    State b1 = basis[1].clone();
    b1.multiplyScalar(df1);
    differential.add(b1);

    double df2 = basis[2].differentiate(fn);
    State b2 = basis[2].clone();
    b2.multiplyScalar(df2);
    differential.add(b2);

    // #//now the differential needs to be converted from a covector to a vector
    // #//using the inverse metric:
    Eigen::Matrix3d metric = this->metricTensor(pos);
    if(metric.determinant() < 0.00001)
    {
        std::cout << "Warning! Metric Tensor Near Singular" << std::endl;
    }

    Eigen::Matrix3d invMetric = metric.inverse();
    differential.vel = invMetric * differential.vel;

    return differential;
}

// This version overloaded to handle distance between vertices in configuration space class
State Geometry::gradient(std::function<double(Eigen::Vector3d,Eigen::Vector3d)> fn, Eigen::Vector3d posi, Eigen::Vector3d pos)
{
    std::vector<State> basis = tangentBasis(pos);
    State differential = State(pos, Eigen::Vector3d(0.,0.,0.));

    // #//add them all up:
    double df0 = basis[0].differentiate(fn, posi);
    State b0 = basis[0].clone();
    b0.multiplyScalar(df0);
    differential.add(b0);

    double df1 = basis[1].differentiate(fn, posi);
    State b1 = basis[1].clone();
    b1.multiplyScalar(df1);
    differential.add(b1);

    double df2 = basis[2].differentiate(fn, posi);
    State b2 = basis[2].clone();
    b2.multiplyScalar(df2);
    differential.add(b2);

    // #//now the differential needs to be converted from a covector to a vector
    // #//using the inverse metric:
    Eigen::Matrix3d metric = this->metricTensor(pos);
    if(metric.determinant() < 0.00001)
    {
        std::cout << "Warning! Metric Tensor Near Singular" << std::endl;
    }

    Eigen::Matrix3d invMetric = metric.inverse();
    differential.vel = invMetric * differential.vel;

    return differential;
}

