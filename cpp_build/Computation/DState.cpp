#include <vector>
#include <iostream>
#include "State.h"
#include "DState.h"


//build a dstate from the input of an object storing velocity data and an object storing acceleration data
DState::DState(std::vector<double> vel, std::vector<double> acc)
{
    this->vel = vel;
    this->acc = acc;
}

//make a copy of a given dstate (not just reference it in memory) and return the copy
DState DState::clone()
{
    return  DState(this->vel,this->acc);
}

//print information about the dstate to the terminal
void DState::print_info()
{
    
    std::cout << "DState Position: [" << this->vel[0] << " , " << this->vel[1] << " , " << this->vel[2] << " ]\n" << "DState Acceleration: [" << this->acc[0] << " , " << this->acc[1] << " , " << this->acc[2] << " ]\n";
}

//add the velocity AND acceleration of a given dstate to the current
void DState::add(DState dstate)
{ 
    // Add the elements of state.vec and this->vel
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        this->vel[i] = this->vel[i] + dstate.vel[i];
        this->acc[i] = this->acc[i] + dstate.acc[i]; 
    }
}

//substract the velocity AND acceleration of a given dstate from the current
void DState::sub(DState dstate)
{
    // Subtract the elements of state.vel from this->vel
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        this->vel[i] = this->vel[i] - dstate.vel[i];
        this->acc[i] = this->acc[i] - dstate.acc[i]; 
    }
}

//scale the velocity AND acceleration of the current state by a factor
void DState::multiplyScalar(double k)
{
    // Multiply the elements of state.vel by k
    for (int i = 0; i < this->vel.size(); i++)  
    {  
        this->vel[i] = this->vel[i] * k; 
        this->acc[i] = this->acc[i] * k;  
    }
}
