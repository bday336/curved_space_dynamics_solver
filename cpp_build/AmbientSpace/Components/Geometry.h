#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <typeinfo>
#include <variant>
#include "../../includes/eigen-3.4.0/Eigen/Dense"
#include "../../Computation/State.h"
#include "../../Computation/DState.h"

// AmbientSpace Class Declaration

// Struct to container helper derivative functions from user input file
struct derivative_funcs {

    // Generate Jacobian for implicit solver methods
    std::function<Eigen::Matrix<double,6,6>(State state)> vertex_jacobian;

    // Pariwise coupling potential derivatives (harmonic) - array of derivatives
    std::function<Eigen::Vector<double,6>(State state1, State state2)> coupling_derviative1_vec;

    // Pariwise coupling potential double derivatives (harmonic) - array of double derivatives
    std::function<Eigen::Matrix<double,6,6>(State state1, State state2)> coupling_derviative2_mat;

    derivative_funcs(){}
    derivative_funcs(
        std::function<Eigen::Matrix<double,6,6>(State state)> vertex_jacobian,
        std::function<Eigen::Vector<double,6>(State state1, State state2)> coupling_derviative1_vec,
        std::function<Eigen::Matrix<double,6,6>(State state1, State state2)> coupling_derviative2_mat){
            this->vertex_jacobian = vertex_jacobian;
            this->coupling_derviative1_vec = coupling_derviative1_vec;
            this->coupling_derviative2_mat = coupling_derviative2_mat;

    }

};

class Geometry
{
    public:

        // Properties
        std::function< Eigen::Matrix3d (Eigen::Vector3d) > metricTensor;
        std::function< Eigen::Vector3d (State) > christoffel;
        std::function< double (Eigen::Vector3d,Eigen::Vector3d)> distance;
        struct derivative_funcs funcDict;

        // Constructor
        Geometry();
        Geometry(
            std::function< Eigen::Matrix3d (Eigen::Vector3d)> metricTensor, 
            std::function< Eigen::Vector3d (State)> christoffel, 
            std::function<double(Eigen::Vector3d,Eigen::Vector3d)> distance
            );
        Geometry(
            std::function< Eigen::Matrix3d (Eigen::Vector3d)> metricTensor, 
            std::function< Eigen::Vector3d (State)> christoffel, 
            std::function<double(Eigen::Vector3d,Eigen::Vector3d)> distance, 
            struct derivative_funcs funcDict
            );

        //Methods

        DState covariantAcceleration(State state);

        double dot(State state1, State state2);

        std::vector<State> tangentBasis(Eigen::Vector3d pos);

        State gradient(std::function<double(Eigen::Vector3d)> fn, Eigen::Vector3d pos);
        // This version overloaded to handle distance between vertices in configuration space class
        State gradient(std::function<double(Eigen::Vector3d,Eigen::Vector3d)> fn, Eigen::Vector3d posi, Eigen::Vector3d pos);

        // methods acting as proxy calls from the 

};


#endif