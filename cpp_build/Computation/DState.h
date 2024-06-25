#ifndef DSTATE_H
#define DSTATE_H

#include <vector>
#include "../includes/eigen-3.4.0/Eigen/Dense"
#include "State.h"

//The class "dState" stores a differential of a
// tangent vector to configuration space
//in practice, this will the configuration space of a single ball system,
//meaning that state.vel will be the ball's velocity and state.acc will be
//covariant derivative with respect to the metric on the ambient space
//so long as they implement the following methods:

class DState
{
    public:

        // Properties
        Eigen::Vector3d vel;
        Eigen::Vector3d acc;

        // std::vector<double> vel;
        // std::vector<double> acc;

        // Constructor
        //build a state from the input of an object storing velocity data and an object storing acceleration data
        DState(Eigen::Vector3d vel, Eigen::Vector3d acc);

        // DState(std::vector<double> vel, std::vector<double> acc);

        // Methods
        //make a copy of a given dstate (not just reference it in memory) and return the copy
        DState clone();

        //print information about the dstate to the terminal
        void print_info();

        //add the velocity AND acceleration of a given state to the current
        void add(DState dstate);

        //subtract the velocity AND acceleration of a given state from the current
        void sub(DState dstate);

        //scale the velocity AND acceleration of the current state by a factor
        void multiplyScalar(double k);


};

#endif