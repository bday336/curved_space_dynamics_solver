#include <vector>
#include "../../includes/eigen-3.4.0/Eigen/Dense"
#include "../../AmbientSpace/AmbientSpace.h"
#include "../State.h"
#include "../DState.h"
#include "../DataList.h"

#include "Solver.h"


//Constructor
Solver::Solver()
{

}
Solver::Solver(AmbientSpace ambientSpace, double stepSize)
{
    this->ambientSpace = ambientSpace;
    this->stepSize = stepSize;
}

//Methods
