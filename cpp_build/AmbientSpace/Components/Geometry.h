#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <typeinfo>
#include <variant>
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

class Geometry
{
    public:

        // Properties
        std::function<std::vector<std::vector<double>>(std::vector<double>)> metricTensor;
        std::function<std::vector<double>(State)> christoffel;
        std::function<double(std::vector<double>,std::vector<double>)> distance;
        GeometryData funcDict;

        // Constructor
        Geometry(
            std::function<std::vector<std::vector<double>>(std::vector<double>)> metricTensor, 
            std::function<std::vector<double>(State)> christoffel, 
            std::function<double(std::vector<double>,std::vector<double>)> distance, 
            GeometryData funcDict
            );

        //Methods

        DState covariantAcceleration(State state);

        double dot(State state1, State state2);

        std::vector<double> tangentBasis(std::vector<double> pos);

        State gradient(std::function<double(std::vector<double>)> fn, std::vector<double> pos);

};


#endif