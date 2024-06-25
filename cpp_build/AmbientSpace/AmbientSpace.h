#ifndef AMBIENT_SPACE_H
#define AMBIENT_SPACE_H

#include <vector>
#include "../includes/eigen-3.4.0/Eigen/Dense"
#include "Components/Geometry.h"
#include "Components/Model.h"
#include "Components/Obstacle.h"
#include "../Computation/State.h"
#include "../Computation/DState.h"

// AmbientSpace Class Declaration

struct ModelData {
    Eigen::Vector3d pos;
    double scaling;
};

class AmbientSpace
{
    public:
        // Properties
        Geometry geometry;
        Model model;
        Obstacle obstacle;


        // Constructor
        AmbientSpace();
        AmbientSpace(Geometry geometry, Model model, Obstacle obstacle);

        //Methods

        DState acceleration(State state);

        double dot(State state1, State state2);

        State gradient(std::function<double(Eigen::Vector3d)> fn, Eigen::Vector3d pos);

        // This version overloaded to handle distance between vertices in configuration space class
        State gradient(std::function<double(Eigen::Vector3d,Eigen::Vector3d)> fn, Eigen::Vector3d posi, Eigen::Vector3d pos);

        ModelData toR3(Eigen::Vector3d pos);

        double distToObstacle(Eigen::Vector3d pos);





};


    // def __init__(self, geometry, model, obstacle):
    //     self.geometry = geometry
    //     self.model    = model
    //     self.obstacle = obstacle

    // # For adding force field
    // def acceleration(self, state):
    //     return self.geometry.covariantAcceleration(state)

    // # Dot product of geometry
    // def dot(self, state1, state2):
    //     return self.geometry.dot(state1,state2)

    // # Gradient with respect to geometry
    // def gradient(self, fn, pos):
    //     return self.geometry.gradient(fn,pos)

    // # Projection from model for visualization
    // def toR3(self, pos):
    //     posR3 = self.model.toR3(pos)
    //     return {pos: posR3, scaling: self.model.relativeScaling(posR3)}

    // # Distance function of geometry
    // def distance(self, pos1, pos2):
    //     return self.geometry.distance(pos1, pos2)

    // # Distance function to the obstacle / bounding box of the simulation
    // def distToObstacle(self, pos):
    //     return self.obstacle.distance(pos)


#endif

