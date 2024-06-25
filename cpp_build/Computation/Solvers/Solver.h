#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include "../../includes/eigen-3.4.0/Eigen/Dense"
#include "../../AmbientSpace/AmbientSpace.h"
#include "../State.h"
#include "../DState.h"
#include "../DataList.h"

class Solver
{
    public:
        //Properties
        AmbientSpace ambientSpace;
        double stepSize;

        //temp variables TO-DO make param argument
        double ks = 1.;
        double x = 1.;
        double m = 1.;

    //Constructor
    Solver();
    Solver(AmbientSpace ambientSpace, double stepSize);

    //Methods
};

#endif