#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include "../../includes/eigen-3.4.0/Eigen/Dense"

// AmbientSpace Class Declaration

class Model
{
    public:

        // Properties
        std::function< Eigen::Vector3d (Eigen::Vector3d) > toR3;
        std::function< double (Eigen::Vector3d) > relativeScaling;

        // Constructor
        Model();
        Model(
            std::function< Eigen::Vector3d (Eigen::Vector3d) > toR3,
            std::function< double (Eigen::Vector3d) > relativeScaling
            );

        //Methods

};


#endif

