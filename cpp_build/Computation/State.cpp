#include <vector>
#include <iostream>
#include "State.h"
#include "DState.h"



//build a state from the input of an object storing position data and an object storing velocity data
State::State(std::vector<double> pos, std::vector<double> vel)
{
    this->pos = pos;
    this->vel = vel;
}

//make a copy of a given state (not just reference it in memory) and return the copy
State State::clone()
{
    return  State(this->pos,this->vel);
}

//print information about the state to the terminal
void State::print_info()
{
    
    std::cout << "State Position: [" << this->pos[0] << " , " << this->pos[1] << " , " << this->pos[2] << " ]\n" << "State Velocity: [" << this->vel[0] << " , " << this->vel[1] << " , " << this->vel[2] << " ]\n";
}

//add the velocity of a given state to the current
void State::add(State state)
{ 
    // Add the elements of state.vec and this->vel
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        this->vel[i] = this->vel[i] + state.vel[i];  
    };
}

//subtract the velocity of a given state from the current
void State::sub(State state)
{
    // Subtract the elements of state.vel from this->vel
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        this->vel[i] = this->vel[i] - state.vel[i];  
    }
}

//scale the velocity of the current state by a factor
void State::multiplyScalar(double k)
{
    // Multiply the elements of state.vel by k
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        this->vel[i] = this->vel[i] * k;  
    }
}

//take the directional derivative of a function fn at pos in direction vel:
double State::differentiate(std::function<double(std::vector<double>)> fn)
{
    double eps = 0.00001;
    std::vector<double> pos1, pos2;

    // Populate vectors pos1 and pos2
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        pos1.push_back(this->pos[i] + this->vel[i] * -eps/2.); 
        pos2.push_back(this->pos[i] + this->vel[i] * eps/2.);
    }

    double dval = fn(pos2) - fn(pos1);
    return  dval/eps;
}

//move a state infintesimally along its tangent direction
void State::flow(double eps)
{
    // Add the elements of this->pos and eps*this->vel
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        this->pos[i] = this->pos[i] + eps * this->vel[i];  
    }
}

//update a state (a tangent vector) by infinitesimally flowing along a differential to the state: a pair dState of a velocity and acceleration
void State::updateBy(DState dstate)
{
    // Add the elements of this->pos and this->vel with dstate.vel and dstate.acc, respectively
    for (int i = 0; i < this->vel.size(); i++)  
    {   
        this->pos[i] = this->pos[i] + dstate.vel[i];
        this->vel[i] = this->vel[i] + dstate.acc[i];  
    }
}