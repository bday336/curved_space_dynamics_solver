#ifndef CONFIGURATION_SPACE_H
#define CONFIGURATION_SPACE_H

#include <vector>
#include "../includes/eigen-3.4.0/Eigen/Dense"
#include "../AmbientSpace/Components/Geometry.h"
#include "../AmbientSpace/Components/Model.h"
#include "../AmbientSpace/Components/Obstacle.h"
#include "../AmbientSpace/AmbientSpace.h"
#include "../Computation/State.h"
#include "../Computation/DState.h"
#include "../Computation/DataList.h"



class ConfigurationSpace
{
    public:
        std::vector<double> masses;
        std::vector<double> radii;
        double N = masses.size();
        AmbientSpace ambientspace;

        //Constructor
        ConfigurationSpace();
        ConfigurationSpace(std::vector<double> masses, std::vector<double> radii, AmbientSpace ambientspace);

        //Methods
        double dot(DataList<State> dataList1, DataList<State> dataList2);

        double norm(DataList<State> dataList);

        DataList<State> normalize(DataList<State> dataList);

        std::vector<int> obstacleCollisions(DataList<State> dataList);

        DataList<State> obstacleGradient(DataList<State> dataList, std::vector<int> indices);

        std::vector<std::vector<int>> ballCollisions(DataList<State> dataList);

        DataList<State> ballGradient(DataList<State> dataList, std::vector<std::vector<int>> indices);

        DataList<State> boundaryNormal(DataList<State> dataList, std::vector<std::vector<int>> ball_collisionIndices, std::vector<int> obstacle_collisionIndices);

        DataList<State> reflectIn(DataList<State> dataList, DataList<State> normal);




};


#endif
