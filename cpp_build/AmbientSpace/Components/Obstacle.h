#ifndef OBSTACLE_H
#define OBSTACLE_H

#include <vector>
#include "../../includes/eigen-3.4.0/Eigen/Dense"

// AmbientSpace Class Declaration

class Obstacle
{
    public:

        // Properties
        std::function< double (Eigen::Vector3d) > distance;
        double size;

        // Constructor
        Obstacle();
        Obstacle(
            std::function< double (Eigen::Vector3d) > distance,
            double size
            );

        //Methods

};


#endif

