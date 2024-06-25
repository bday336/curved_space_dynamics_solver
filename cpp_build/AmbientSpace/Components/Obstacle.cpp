#include <vector>
#include "Obstacle.h"
#include "../../includes/eigen-3.4.0/Eigen/Dense"

// AmbientSpace Class Declaration

// Constructor
Obstacle::Obstacle()
{
    
}

Obstacle::Obstacle(
    std::function< double (Eigen::Vector3d) > distance,
    double size
    )
    {
        this->distance = distance;
        this->size = size;
    }

//Methods

