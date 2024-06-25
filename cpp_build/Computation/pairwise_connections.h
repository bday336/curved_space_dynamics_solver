#ifndef PAIRWISE_CONNECTIONS_H
#define PAIRWISE_CONNECTIONS_H

#include <vector>

// Used for pairwise coupling between vertices (harmonic potential)
struct Pairwise_Connection
{
    // Properties
    // Eigen::Vector2i index_pair;

    std::vector<int> index_pair;
};

// Used for pairwise rigid constraint between vertices
struct Pairwise_Rigid_Connection
{
    // Properties
//    Eigen::Vector2i index_pair;
    std::vector<int> index_pair;
    double lagrange_multipler;
};

#endif