#ifndef STATE_H
#define STATE_H

#include <vector>
#include "DState.h"

class State
{
    public:

        // Properties
        std::vector<double> pos;
        std::vector<double> vel;

        // Constructor
        //build a state from the input of an object storing position data and an object storing velocity data
        State(std::vector<double> pos, std::vector<double> vel);

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
        double differentiate(std::function<double(std::vector<double>)> fn);

        //move a state infintesimally along its tangent direction
        void flow(double eps);

        //update a state (a tangent vector) by infinitesimally flowing along a differential to the state: a pair dState of a velocity and acceleration
        void updateBy(DState dstate);

};

#endif