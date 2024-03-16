#include <vector>
#include <typeinfo>
#include <variant>
#include "AmbientSpace.h"
#include "../Computation/State.h"
#include "../Computation/DState.h"

// AmbientSpace Class Declaration

struct ModelData {
    std::vector<double> pos;
    double scaling;
};

AmbientSpace::AmbientSpace(Geometry geometry, Model model, Obstacle obstacle)
{

}

        //Methods

        DState acceleration(State state);

        double dot(State state1, State state2);

        State gradient(std::function<double(std::vector<double>)> fn, std::vector<double> pos)

        ModelData toR3(std::vector<double> pos);

        double distance(std::vector<double> pos1, std::vector<double> pos2);

        double distToObstacle(std::vector<double> pos);





}


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


