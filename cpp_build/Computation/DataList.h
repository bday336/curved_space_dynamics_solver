#ifndef DATALIST_H
#define DATALIST_H

#include <vector>
#include <typeinfo>
#include <variant>
#include "State.h"
#include "DState.h"
#include "pairwise_connections.h"

template <typename T>
class DataList
{
    public:

        // Properties
        std::vector<T> data;
        std::vector<Pairwise_Connection> connectivity;
        std::vector<Pairwise_Rigid_Connection> rig_connectivity;

        // Constructor
        // Populate datalist with states and mesh connectivity data
        DataList(std::vector<T> data, std::vector<Pairwise_Connection> connectivity, std::vector<Pairwise_Rigid_Connection> rig_connectivity);

        // Methods
        DataList clone();

        DataList combine(DataList datalist);

        std::vector<double> toArray();

        //print information about the state to the terminal
        void print_info();

        //add the velocity of a given state to the current
        void add(DataList datalist);

        //subtract the velocity of a given state from the current
        void sub(DataList datalist);

        //scale the velocity of the current state by a factor
        void multiplyScalar(double k);

        //move a state infintesimally along its tangent direction
        void flow(double eps);

        //update a state (a tangent vector) by infinitesimally flowing along a differential to the state: a pair dState of a velocity and acceleration
        void updateBy(DataList<DState> datalist);
};


#endif