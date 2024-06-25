#include <vector>
#include "Model.h"
#include "../../includes/eigen-3.4.0/Eigen/Dense"

// AmbientSpace Class Declaration

// Constructor
Model::Model()
{
    
}

Model::Model(
    std::function< Eigen::Vector3d (Eigen::Vector3d) > toR3,
    std::function< double (Eigen::Vector3d) > relativeScaling
    )
    {
        this->toR3 = toR3;
        this->relativeScaling = relativeScaling;
    }

