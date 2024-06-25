#ifndef STATE_H
#define STATE_H

#include <vector>
#include "../includes/eigen-3.4.0/Eigen/Dense"

#include "DState.h"

class State
{
    public:

        // Properties
        Eigen::Vector3d pos;
        Eigen::Vector3d vel;

        // std::vector<double> pos;
        // std::vector<double> vel;

        // Constructor
        //build a state from the input of an object storing position data and an object storing velocity data
        State(Eigen::Vector3d pos, Eigen::Vector3d vel);

        // State(std::vector<double> pos, std::vector<double> vel);

        // Methods
        //make a copy of a given state (not just reference it in memory) and return the copy
        State clone();

        //print information about the state to the terminal
        void print_info();

        //add the velocity of a given state to the current
        void add(State state);

        //subtract the velocity of a given state from the current
        void sub(State state);

        //scale the velocity of the current state by a factor
        void multiplyScalar(double k);

        //take the directional derivative of a function at pos in direction vel:
        double differentiate(std::function<double(Eigen::Vector3d)> fn);
        // This version overloaded to handle distance between vertices in configuration space class
        double differentiate(std::function<double(Eigen::Vector3d,Eigen::Vector3d)> fn, Eigen::Vector3d posi);

        // double differentiate(std::function<double(std::vector<double>)> fn);

        //move a state infintesimally along its tangent direction
        void flow(double eps);

        //update a state (a tangent vector) by infinitesimally flowing along a differential to the state: a pair dState of a velocity and acceleration
        void updateBy(DState dstate);

};

#endif