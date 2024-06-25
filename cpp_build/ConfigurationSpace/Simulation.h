#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "../includes/eigen-3.4.0/Eigen/Dense"
#include "../AmbientSpace/Components/Geometry.h"
#include "../AmbientSpace/Components/Model.h"
#include "../AmbientSpace/Components/Obstacle.h"
#include "../AmbientSpace/AmbientSpace.h"
#include "../Computation/State.h"
#include "../Computation/DState.h"
#include "../Computation/DataList.h"

#include "ConfigurationSpace.h"



class Simulation
{
    public:
        AmbientSpace ambientSpace;
        DataList<State> dataList;
        ConfigurationSpace configurationSpace;

        // Solver solver_method;
        double stepSize;

        //to set when intersecting
        std::vector<std::vector<int>> ball_collisions;
        std::vector<int >obstacle_collisions;

        //Container for all simulation data
        std::vector<DataList<State>> data_container;


        //build an integrator
        //get the function which takes the derivative of each element of a stateList:
        //using ambientSpace.acceleration will allow us to use external potentials without changing code

        // if self.solver_method == "RungeKutta":
        //     self.integrator = RungeKutta(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Gauss1":
        //     self.integrator = Gauss1(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Gauss1test":
        //     self.integrator = Gauss1test(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Gauss2":
        //     self.integrator = Gauss2(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Gauss2test":
        //     self.integrator = Gauss2test(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Gauss3":
        //     self.integrator = Gauss3(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Radau2":
        //     self.integrator = Radau2(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Radau2test":
        //     self.integrator = Radau2test(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "RigidRadau2":
        //     self.integrator = RigidRadau2(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "Radau3":
        //     self.integrator = Radau3(self.ambientSpace,self.stepSize)
        // elif self.solver_method == "RigidRadau3":
        //     self.integrator = RigidRadau3(self.ambientSpace,self.stepSize)

        //Constructor
        Simulation();
        Simulation(AmbientSpace ambientSpace,DataList<State> dataList,ConfigurationSpace configurationSpace,double stepSize,Solver solver_method);

        //Methods
        bool detectCollision();

        void smoothDynamics();

        void collisionDynamics();

        void step();

};

#endif
