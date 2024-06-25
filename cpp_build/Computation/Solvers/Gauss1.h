#ifndef GAUSS1_H
#define GAUSS1_H

#include <vector>
#include "../../includes/eigen-3.4.0/Eigen/Dense"
#include "../../AmbientSpace/AmbientSpace.h"
#include "../State.h"
#include "../DState.h"
#include "../DataList.h"

#include "Solver.h"

class Gauss1 : public Solver
{
    public:

    //Constructor
    Gauss1();
    Gauss1(AmbientSpace ambientSpace, double stepSize);

    //Methods
    DataList<State> arrayToDataList(Eigen::Vector3d arr, DataList<State> dataList);

    DataList<DState> dynfunc(DataList<State> dataList);

    Eigen::MatrixXd dynjac(DataList<State> dataList);

    Eigen::VectorXd difffunc(DataList<State> dataList, DataList<DState> k1);

    DataList<State> step(DataList<State> dataList, double tol = 1e-15, int imax = 100);

};

#endif