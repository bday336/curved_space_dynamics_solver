#include <vector>
#include "../includes/eigen-3.4.0/Eigen/Dense"
#include "Components/Geometry.h"
#include "Components/Model.h"
#include "Components/Obstacle.h"
#include "AmbientSpace.h"
#include "../Computation/State.h"
#include "../Computation/DState.h"

// AmbientSpace Class Declaration

AmbientSpace::AmbientSpace()
{

}

AmbientSpace::AmbientSpace(Geometry geometry, Model model, Obstacle obstacle)
{
    this->geometry = geometry;
    this->model = model;
    this->obstacle = obstacle;
}

//Methods

DState AmbientSpace::acceleration(State state)
{
    return this->geometry.covariantAcceleration(state);
}

double AmbientSpace::dot(State state1, State state2)
{
    return this->geometry.dot(state1,state2);
}

State AmbientSpace::gradient(std::function<double(Eigen::Vector3d)> fn, Eigen::Vector3d pos)
{
    return this->geometry.gradient(fn,pos);
}

// This version overloaded to handle distance between vertices in configuration space class
State AmbientSpace::gradient(std::function<double(Eigen::Vector3d,Eigen::Vector3d)> fn, Eigen::Vector3d posi, Eigen::Vector3d pos)
{
    return this->geometry.gradient(fn,posi,pos);
}

ModelData AmbientSpace::toR3(Eigen::Vector3d pos)
{
    Eigen::Vector3d posR3 = this->model.toR3(pos);
    ModelData mdat;
    mdat.pos = posR3;
    mdat.scaling = this->model.relativeScaling(posR3);
    return mdat;
}

double AmbientSpace::distToObstacle(Eigen::Vector3d pos)
{
    return this->obstacle.distance(pos);
}


