#ifndef H3_HYP_HPP
#define H3_HYP_HPP

#include <vector>
#include <iostream>
#include <ginac/ginac.h>
#include "../../includes/eigen-3.4.0/Eigen/Dense"
#include "../Components/Geometry.h"
#include "../Components/Model.h"
#include "../Components/Obstacle.h"
#include "../AmbientSpace.h"
#include "../../Computation/State.h"
#include "../../Computation/DState.h"

const double pi = 3.14159265358979323846;

// Define the variables used in symbolic input expressions

GiNaC::symbol a("a"), b("b"), g("g");
GiNaC::symbol ad("ad"), bd("bd"), gd("gd");
GiNaC::symbol a1("a1"), b1("b1"), g1("g1");
GiNaC::symbol a2("a2"), b2("b2"), g2("g2");

// -------------------------------------------------------------
// Geometry Information
// -------------------------------------------------------------

// User provided information

// Coordinates

GiNaC::lst coords = {
    GiNaC::sinh(a)*GiNaC::sin(b)*GiNaC::cos(g),
    GiNaC::sinh(a)*GiNaC::sin(b)*GiNaC::sin(g),
    GiNaC::sinh(a)*GiNaC::cos(b),
    GiNaC::cosh(a)
};

// Metric

GiNaC::matrix metric = {
    {1.,0.,0.},
    {0.,GiNaC::pow(GiNaC::sinh(a),2),0.},
    {0.,0.,GiNaC::pow(GiNaC::sinh(a),2)*GiNaC::pow(GiNaC::sin(b),2)}
};

// Distance

GiNaC::ex distance = GiNaC::acosh(
    -GiNaC::sinh(a1)*GiNaC::sin(b1)*GiNaC::cos(g1)*GiNaC::sinh(a2)*GiNaC::sin(b2)*GiNaC::cos(g2) -
     GiNaC::sinh(a1)*GiNaC::sin(b1)*GiNaC::sin(g1)*GiNaC::sinh(a2)*GiNaC::sin(b2)*GiNaC::sin(g2) -
     GiNaC::sinh(a1)*GiNaC::cos(b1)*GiNaC::sinh(a2)*GiNaC::cos(b2) +
     GiNaC::cosh(a1)*GiNaC::cosh(a2)
    );

// Pairwise coupling potential (harmonic)

GiNaC::symbol m("m"), k("k"), lo("lo");

GiNaC::ex pairwise_potential = .5*k*GiNaC::pow((distance - lo),2);

// Derived expressions for numerics using GiNaC

// Metric input function
Eigen::Matrix3d geoMetric(Eigen::Vector3d pos){
    Eigen::Matrix3d mat;
    GiNaC::exmap gmap;
    gmap[a] = pos[0];
    gmap[b] = pos[1];
    gmap[g] = pos[2];
    for (int i = 0; i < mat.rows(); i++)
    {
       for (int j = 0; j < mat.cols(); j++)
       {
        mat(i,j) = GiNaC::ex_to<GiNaC::numeric>(metric(i,j).subs(gmap)).to_double();
       }
       
    }
    return mat;
};

GiNaC::matrix inverse_metric = metric.inverse();

GiNaC::symbol isym("isym"), jsym("jsym"), ksym("ksym");

GiNaC::ex christFunc(int i, int j, int k){
    switch(i){
        case 0:
            isym = a;
            break;
        case 1:
            isym = b;
            break;
        case 2:
            isym = g;
            break;
    }
    switch(j){
        case 0:
            jsym = a;
            break;
        case 1:
            jsym = b;
            break;
        case 2:
            jsym = g;
            break;
    }
    switch(k){
        case 0:
            ksym = a;
            break;
        case 1:
            ksym = b;
            break;
        case 2:
            ksym = g;
            break;
    }
    GiNaC::ex christ = .5*inverse_metric(i,0)*(metric(0,j).diff(ksym)+metric(0,k).diff(jsym)-metric(j,k).diff(isym)) +
                       .5*inverse_metric(i,1)*(metric(1,j).diff(ksym)+metric(1,k).diff(jsym)-metric(j,k).diff(isym)) +
                       .5*inverse_metric(i,2)*(metric(2,j).diff(ksym)+metric(2,k).diff(jsym)-metric(j,k).diff(isym));
    return christ;
}

// Negative since equals second derivative idd = -G^i_JK * jd * kd
GiNaC::lst christoffel = {
    -1.*(ad*ad*christFunc(0,0,0)+ad*bd*christFunc(0,0,1)+ad*gd*christFunc(0,0,2)+
         bd*ad*christFunc(0,1,0)+bd*bd*christFunc(0,1,1)+bd*gd*christFunc(0,1,2)+
         gd*ad*christFunc(0,2,0)+gd*bd*christFunc(0,2,1)+gd*gd*christFunc(0,2,2)),

    -1.*(ad*ad*christFunc(1,0,0)+ad*bd*christFunc(1,0,1)+ad*gd*christFunc(1,0,2)+
         bd*ad*christFunc(1,1,0)+bd*bd*christFunc(1,1,1)+bd*gd*christFunc(1,1,2)+
         gd*ad*christFunc(1,2,0)+gd*bd*christFunc(1,2,1)+gd*gd*christFunc(1,2,2)),

    -1.*(ad*ad*christFunc(2,0,0)+ad*bd*christFunc(2,0,1)+ad*gd*christFunc(2,0,2)+
         bd*ad*christFunc(2,1,0)+bd*bd*christFunc(2,1,1)+bd*gd*christFunc(2,1,2)+
         gd*ad*christFunc(2,2,0)+gd*bd*christFunc(2,2,1)+gd*gd*christFunc(2,2,2))
    };

// Christoffel input function
Eigen::Vector3d geoChristoffel(State state){
    Eigen::Vector3d vec;
    GiNaC::exmap gmap;
    gmap[a] = state.pos[0];
    gmap[b] = state.pos[1];
    gmap[g] = state.pos[2];

    gmap[ad] = state.vel[0];
    gmap[bd] = state.vel[1];
    gmap[gd] = state.vel[2];
    for (int i = 0; i < vec.size(); i++)
    {
        vec(i) = GiNaC::ex_to<GiNaC::numeric>(christoffel[i].subs(gmap)).to_double();       
    }
    return vec;
};

// For Euler-Lagrange equations
GiNaC::lst el_d1_pairwise_potential = {
    pairwise_potential.diff(a1)/(m*metric(0,0).subs(GiNaC::lst {a == a1, b == b1, g == g1})),
    pairwise_potential.diff(b1)/(m*metric(1,1).subs(GiNaC::lst {a == a1, b == b1, g == g1})),
    pairwise_potential.diff(g1)/(m*metric(2,2).subs(GiNaC::lst {a == a1, b == b1, g == g1})),

    pairwise_potential.diff(a2)/(m*metric(0,0).subs(GiNaC::lst {a == a2, b == b2, g == g2})),
    pairwise_potential.diff(b2)/(m*metric(1,1).subs(GiNaC::lst {a == a2, b == b2, g == g2})),
    pairwise_potential.diff(g2)/(m*metric(2,2).subs(GiNaC::lst {a == a2, b == b2, g == g2}))
};

// For jacobian for implicit solver
GiNaC::matrix el_d2_pairwise_potential = {
    {
        el_d1_pairwise_potential[0].diff(a1),
        el_d1_pairwise_potential[1].diff(a1),
        el_d1_pairwise_potential[2].diff(a1),
        el_d1_pairwise_potential[3].diff(a1),
        el_d1_pairwise_potential[4].diff(a1),
        el_d1_pairwise_potential[5].diff(a1)
    },
    {
        el_d1_pairwise_potential[0].diff(b1),
        el_d1_pairwise_potential[1].diff(b1),
        el_d1_pairwise_potential[2].diff(b1),
        el_d1_pairwise_potential[3].diff(b1),
        el_d1_pairwise_potential[4].diff(b1),
        el_d1_pairwise_potential[5].diff(b1)
    },
    {
        el_d1_pairwise_potential[0].diff(g1),
        el_d1_pairwise_potential[1].diff(g1),
        el_d1_pairwise_potential[2].diff(g1),
        el_d1_pairwise_potential[3].diff(g1),
        el_d1_pairwise_potential[4].diff(g1),
        el_d1_pairwise_potential[5].diff(g1)
    },
    {
        el_d1_pairwise_potential[0].diff(a2),
        el_d1_pairwise_potential[1].diff(a2),
        el_d1_pairwise_potential[2].diff(a2),
        el_d1_pairwise_potential[3].diff(a2),
        el_d1_pairwise_potential[4].diff(a2),
        el_d1_pairwise_potential[5].diff(a2)
    },
    {
        el_d1_pairwise_potential[0].diff(b2),
        el_d1_pairwise_potential[1].diff(b2),
        el_d1_pairwise_potential[2].diff(b2),
        el_d1_pairwise_potential[3].diff(b2),
        el_d1_pairwise_potential[4].diff(b2),
        el_d1_pairwise_potential[5].diff(b2)
    },
    {
        el_d1_pairwise_potential[0].diff(g2),
        el_d1_pairwise_potential[1].diff(g2),
        el_d1_pairwise_potential[2].diff(g2),
        el_d1_pairwise_potential[3].diff(g2),
        el_d1_pairwise_potential[4].diff(g2),
        el_d1_pairwise_potential[5].diff(g2)
    }
};

GiNaC::matrix christMat = {
    {0.,0.,0., 1.,0.,0.},
    {0.,0.,0., 0.,1.,0.},
    {0.,0.,0., 0.,0.,1.},
    {christoffel[0].diff(a),christoffel[0].diff(b),christoffel[0].diff(g), christoffel[0].diff(ad),christoffel[0].diff(bd),christoffel[0].diff(gd)},
    {christoffel[1].diff(a),christoffel[1].diff(b),christoffel[1].diff(g), christoffel[1].diff(ad),christoffel[1].diff(bd),christoffel[1].diff(gd)},
    {christoffel[2].diff(a),christoffel[2].diff(b),christoffel[2].diff(g), christoffel[2].diff(ad),christoffel[2].diff(bd),christoffel[2].diff(gd)}
};

// Distance input function
double geoDistance(Eigen::Vector3d pos1, Eigen::Vector3d pos2){
    GiNaC::exmap gmap;
    gmap[a1] = pos1[0];
    gmap[b1] = pos1[1];
    gmap[g1] = pos1[2];
    gmap[a2] = pos2[0];
    gmap[b2] = pos2[1];
    gmap[g2] = pos2[2];

    return GiNaC::ex_to<GiNaC::numeric>(distance.subs(gmap)).to_double();
};


// Generate Jacobian for implicit solver methods
Eigen::Matrix<double,6,6> vertex_jacobian(State state){
    Eigen::Matrix<double,6,6> mat;
    GiNaC::exmap gmap;
    gmap[a] = state.pos[0];
    gmap[b] = state.pos[1];
    gmap[g] = state.pos[2];

    gmap[ad] = state.vel[0];
    gmap[bd] = state.vel[1];
    gmap[gd] = state.vel[2];
    for (int i = 0; i < mat.rows(); i++)
    {
        for (int j = 0; j < mat.cols(); j++)
        {
            mat(i,j) = GiNaC::ex_to<GiNaC::numeric>(christMat(i,j).subs(gmap)).to_double();
        }
    }
    return mat;
};

// Pariwise coupling potential derivatives (harmonic) - array of derivatives
Eigen::Vector<double,6> coupling_derviative1_vec(State state1, State state2){
    Eigen::Vector<double,6> vec;
    GiNaC::exmap gmap;
    gmap[a1] = state1.pos[0];
    gmap[b1] = state1.pos[1];
    gmap[g1] = state1.pos[2];

    gmap[a2] = state2.pos[0];
    gmap[b2] = state2.pos[1];
    gmap[g2] = state2.pos[2];

    for (int i = 0; i < vec.size(); i++)
    {
        vec(i) = GiNaC::ex_to<GiNaC::numeric>(el_d1_pairwise_potential[i].subs(gmap)).to_double();       
    }
    return vec;
};

// Pariwise coupling potential double derivatives (harmonic) - array of double derivatives
Eigen::Matrix<double,6,6> coupling_derviative2_mat(State state1, State state2){
    Eigen::Matrix<double,6,6> mat;
    GiNaC::exmap gmap;
    gmap[a1] = state1.pos[0];
    gmap[b1] = state1.pos[1];
    gmap[g1] = state1.pos[2];

    gmap[a2] = state2.pos[0];
    gmap[b2] = state2.pos[1];
    gmap[g2] = state2.pos[2];

    for (int i = 0; i < mat.rows(); i++)
    {
        for (int j = 0; j < mat.cols(); j++)
        {
            mat(i,j) = GiNaC::ex_to<GiNaC::numeric>(el_d2_pairwise_potential(i,j).subs(gmap)).to_double();
        }
    }
    return mat;
};

// Helper functions for solver based on the geometry and system data
derivative_funcs geo_dfuncs(
    vertex_jacobian,
    coupling_derviative1_vec,
    coupling_derviative2_mat
);

Geometry geometry = Geometry(
    geoMetric,
    geoChristoffel,
    geoDistance,
    geo_dfuncs
    );

// -------------------------------------------------------------
// Model Information
// -------------------------------------------------------------

Eigen::Vector3d toPoincareBall(Eigen::Vector3d pos){
    Eigen::Vector3d vec;
    GiNaC::exmap gmap;
    gmap[a] = pos[0];
    gmap[b] = pos[1];
    gmap[g] = pos[2];
    for (int i = 0; i < vec.size(); i++)
    {
        vec(i) = GiNaC::ex_to<GiNaC::numeric>(coords[i].subs(gmap)).to_double() / (1. + GiNaC::ex_to<GiNaC::numeric>(coords[3].subs(gmap)).to_double());       
    }
    return vec;
};

double pbScaling(Eigen::Vector3d pos){
    double modelScale = 6.;
    double scale = modelScale*modelScale - pos.squaredNorm();
    return 4. * scale / (modelScale*modelScale);
};

Model model = Model(
    toPoincareBall,
    pbScaling
    );

// -------------------------------------------------------------
// Obstacle/Bounding Ball Information
// -------------------------------------------------------------

// Sphere Bounding Box

// Default bounding box size (radius)
double obstacleSize = 2;

double distToSphere(Eigen::Vector3d pos){
    //center point of H3 in coordinates:
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    //distance from center to position
    double dist = geoDistance(pos,center);
    //how far is this from the boundary sphere of radius 5?
    return obstacleSize - dist;
};


Obstacle sphereObstacle = Obstacle(
    distToSphere,
    obstacleSize
);

//package stuff up for export
AmbientSpace hyperbolic = AmbientSpace(geometry, model, sphereObstacle);

#endif
