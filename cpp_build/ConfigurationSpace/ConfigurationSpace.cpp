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


//Constructor
ConfigurationSpace::ConfigurationSpace()
{

}

ConfigurationSpace::ConfigurationSpace(std::vector<double> masses, std::vector<double> radii, AmbientSpace ambientspace)
{
    this->masses = masses;
    this->radii = radii;
    this->ambientspace = ambientspace;
}

//Methods
double ConfigurationSpace::dot(DataList<State> dataList1, DataList<State> dataList2)
{
    double dot = 0;
    for (int i = 0; i < this->N; i++)
    {
        dot = dot + .5 * this->masses[i] * this->ambientspace.dot(dataList1.data[i],dataList2.data[i]);
    }
    return dot;
    
}

double ConfigurationSpace::norm(DataList<State> dataList)
{
    double val = std::sqrt(this->dot(dataList,dataList));
    return val;
}

DataList<State> ConfigurationSpace::normalize(DataList<State> dataList)
{
    double len = this->norm(dataList);
    DataList<State> res = dataList;
    res.multiplyScalar(1./len);
    return res;
}

std::vector<int> ConfigurationSpace::obstacleCollisions(DataList<State> dataList)
{
    std::vector<int> indices;
    for (int i = 0; i < this->N; i++)
    {
        Eigen::Vector3d posi = dataList.data[i].pos;
        double disti = this->ambientspace.distToObstacle(posi);
        double radi = this->radii[i];
        if (disti < radi)
        {
            //the balls is intersecting the boundary:
            //but, see if it is heading outward or inward
            State newState = dataList.data[i].clone();
            newState.flow(.001);
            Eigen::Vector3d newPos = newState.pos;
            double newDist = this->ambientspace.distToObstacle(newPos);
            //if this new distance is less, it's an intersection with inadmissable tangent
            if (newDist < disti)
            {
                indices.push_back(i);
            }
        }
    }

    return indices;
    
    
}

DataList<State> ConfigurationSpace::obstacleGradient(DataList<State> dataList, std::vector<int> indices)
{
    //make a new state with same positions but zeroed out velocity:
    DataList<State> grad = dataList.clone();
    Eigen::Vector3d vec(0.,0.,0.);
    for (int i = 0; i < this->N; i++)
    {
        grad.data[i].vel = vec;
    }

    if (indices.size() != 0)
    {
        //replace the velocity with the gradient in the correct index slots
        for (int i = 0; i < indices.size(); i++)
        {
            int ival = indices[i];
            Eigen::Vector3d posi = dataList.data[ival].pos;

            //with respect to the metric g on the ambient space X
            State geomGradi = this->ambientspace.gradient(this->ambientspace.obstacle.distance, posi);

            //the kinetic energy metric is 1/2*m*g, so the inverse metric tensor is scaled by 2/m:
            geomGradi.multiplyScalar(2./this->masses[ival]);

            //replace this in the gradient list:
            grad.data[i] = geomGradi;
        }
        
    }

    return grad;
}

std::vector<std::vector<int>> ConfigurationSpace::ballCollisions(DataList<State> dataList)
{
    std::vector<std::vector<int>> indices;
    for (int i = 0; i < this->N; i++)
    {
        for (int j = 0; j < this->N; j++)
        {
            if (i != j)
            {
                double distij = this->ambientspace.geometry.distance(dataList.data[i].pos,dataList.data[j].pos);
                double radij = this->radii[i] + this->radii[j];

                if (distij < radij)
                {
                    //the balls are intersecting: but are they approaching or moving apart?
                    State newStatei = dataList.data[i].clone();
                    newStatei.flow(.001);
                    Eigen::Vector3d newPosi = newStatei.pos;

                    State newStatej = dataList.data[j].clone();
                    newStatej.flow(.001);
                    Eigen::Vector3d newPosj = newStatej.pos;

                    double newDist = this->ambientspace.geometry.distance(newPosi,newPosj);
                    //if this new distance is less, it's an intersection with inadmissable tangent
                    if (newDist < distij)
                    {
                        std::vector<int> idx = {i,j};
                        indices.push_back(idx);
                    }
                }
            }
        }
    }
    return indices;
}

DataList<State> ConfigurationSpace::ballGradient(DataList<State> dataList, std::vector<std::vector<int>> indices)
{
    //make a new state with same positions but zeroed out velocity:
    DataList<State> grad = dataList.clone();
    Eigen::Vector3d vec(0.,0.,0.);
    for (int i = 0; i < this->N; i++)
    {
        grad.data[i].vel = vec;
    }

    //replace the velocity with the gradient in the correct index slots
    if (indices.size() != 0)
    {
        for (int i = 0; i < indices.size(); i++)
        {
            std::vector<int> ij = indices[i];
            int idxi = ij[0];
            int idxj = ij[1];

            Eigen::Vector3d posi = dataList.data[idxi].pos;
            Eigen::Vector3d posj = dataList.data[idxj].pos;

            //the gradient of this function, evaluated at position j
            State gradjdisti = this->ambientspace.gradient(
                this->ambientspace.geometry.distance, posi, posj
            );

            //the kinetic energy metric for the jth particle is the Riemannian metric g,
            //scaled by 1/2*m: thus the gradient is the g-gradient scaled by 2/m
            //replace the gradient at j with this:
            grad.data[idxj] = gradjdisti.clone();
            grad.data[idxj].multiplyScalar(2./this->masses[idxj]);

            //the gradient of this function, evaluated at position j
            State gradidistj = this->ambientspace.gradient(
                this->ambientspace.geometry.distance, posj, posi
            );

            //the kinetic energy metric for the jth particle is the Riemannian metric g,
            //scaled by 1/2*m: thus the gradient is the g-gradient scaled by 2/m
            //replace the gradient at i with this:
            grad.data[idxi] = gradidistj.clone();
            grad.data[idxi].multiplyScalar(2./this->masses[idxi]);
        }
    }
    return grad;
    
}

DataList<State> ConfigurationSpace::boundaryNormal(DataList<State> dataList, std::vector<std::vector<int>> ball_collisionIndices, std::vector<int> obstacle_collisionIndices)
{
    DataList<State> grad1 = this->obstacleGradient(dataList, obstacle_collisionIndices);
    DataList<State> grad2 = this->ballGradient(dataList, ball_collisionIndices);

    //Add together
    grad1.add(grad2);

    return this->normalize(grad1);
}

//reflect a state in a normal vector
//dataList is the current tangent vector to configuration space (dataList of all the balls)
//normal is the normal vector to the boundary of configuration space
DataList<State> ConfigurationSpace::reflectIn(DataList<State> dataList, DataList<State> normal)
{
    double dot = this->dot(dataList,normal);
    double norm2 = this->dot(normal,normal);

    double coef = 2.*dot/norm2;

    dataList.sub(normal);
    dataList.multiplyScalar(coef);

    return dataList;
}

