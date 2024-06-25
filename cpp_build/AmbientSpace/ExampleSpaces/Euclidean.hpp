#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H

#include <vector>
#include <cmath>
#include "../../includes/eigen-3.4.0/Eigen/Dense"

#include "../Components/Geometry.h"
#include "../Components/Model.h"
#include "../Components/Obstacle.h"
#include "../AmbientSpace.h"
#include "../../Computation/State.h"

// -------------------------------------------------------------
// 3-Dimensional Euclidean Space Information
// -------------------------------------------------------------

// Eigen::Matrix<double, 6, 6> euc_gen_jacobian(State state){
//     Eigen::Matrix<double, 6, 6> mat;
//     return mat;
// };

// GeometryData eucGeoData;
// eucGeoData.vertex_jacobian = euc_gen_jacobian;

// def hyp_gen_jacobian(state):

//     pos = state.pos.copy()
//     vel = state.vel.copy()

//     alpha = pos[0]
//     beta  = pos[1]
//     gamma = pos[2]

//     dalpha = vel[0]
//     dbeta  = vel[1]
//     dgamma = vel[2]

//     # Below the labeling corresponds with:
//     # daddalpha  -> derivative of ddalpha (see hypChristoffel()) with respect to alpha
//     # ddaddalpha -> derivative of ddalpha (see hypChristoffel()) with respect to dalpha
//     #---------- 

//     daddalpha = cosh(2.*alpha) * (dbeta**2. + sin(beta)**2.*dgamma**2.)
//     dbddalpha = cos(beta)*sin(beta)*sinh(2.*alpha)*dgamma**2.
//     dgddalpha = 0.

//     ddaddalpha = 0.
//     ddbddalpha = sinh(2.*alpha)*dbeta
//     ddgddalpha = sin(beta)**2.*sinh(2.*alpha)*dgamma

//     #----------

//     daddbeta = 2./sinh(alpha)**2.*(dalpha*dbeta)
//     dbddbeta = cos(2.*beta)*dgamma**2.
//     dgddbeta = 0.

//     ddaddbeta = -2./tanh(alpha)*dbeta
//     ddbddbeta = -2./tanh(alpha)*dalpha
//     ddgddbeta = sin(2.*beta)*dgamma

//     #------------

//     daddgamma = 2./sinh(alpha)**2.*(dalpha*dgamma)
//     dbddgamma = 2./sin(beta)**2.*(dbeta*dgamma)
//     dgddgamma = 0.

//     ddaddgamma = -2./tanh(alpha)*dgamma
//     ddbddgamma = -2./tan(beta)*dgamma
//     ddgddgamma = -2.*(dalpha/tanh(alpha) + dbeta/tan(beta))


//     return np.array([
//         [0., 0., 0.,  1., 0., 0.],
//         [0., 0., 0.,  0., 1., 0.],
//         [0., 0., 0.,  0., 0., 1.],
//         [daddalpha,dbddalpha,dgddalpha, ddaddalpha,ddbddalpha,ddgddalpha],
//         [daddbeta ,dbddbeta ,dgddbeta , ddaddbeta ,ddbddbeta ,ddgddbeta ],
//         [daddgamma,dbddgamma,dgddgamma, ddaddgamma,ddbddgamma,ddgddgamma]
//     ])

// ##### Distance Functions for coupling potentials and length constraints
// def D12(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return np.cosh(a1)*np.cosh(a2) - np.sinh(a1)*np.cos(b1)*np.sinh(a2)*np.cos(b2) - np.sinh(a1)*np.sin(b1)*np.sinh(a2)*np.sin(b2)*np.cos(g1 - g2)

// # First Derivatives of distance function
// def dtD12(dstate1, dstate2, dterms):
//     da1,db1,dg1 = dstate1.vel.copy()
//     da2,db2,dg2 = dstate2.vel.copy()
//     return (np.array([da1,db1,dg1,da2,db2,dg2]) @ np.array(dterms))

// def da1D12(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return np.sinh(a1)*np.cosh(a2) - np.cosh(a1)*np.cos(b1)*np.sinh(a2)*np.cos(b2) - np.cosh(a1)*np.sin(b1)*np.sinh(a2)*np.sin(b2)*np.cos(g1 - g2)

// def db1D12(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return np.sinh(a1)*np.sin(b1)*np.sinh(a2)*np.cos(b2) - np.sinh(a1)*np.cos(b1)*np.sinh(a2)*np.sin(b2)*np.cos(g1 - g2) 

// def dg1D12(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return np.sinh(a1)*np.sin(b1)*np.sinh(a2)*np.sin(b2)*np.sin(g1 - g2)
// # For the remaining three functions use:
// # da2D12 = da1D12(state2, state1)
// # db2D12 = db1D12(state2, state1)
// # dg2D12 = dg1D12(state2, state1)

// # Second Derivatives of distance function (Only needed for implicit methods for jacobian construction)
// def ddtD12(dstate1, dstate2, dterms, ddterms):
//     da1,db1,dg1 = dstate1.vel.copy()
//     da2,db2,dg2 = dstate2.vel.copy()
//     dda1,ddb1,ddg1 = dstate1.acc.copy()
//     dda2,ddb2,ddg2 = dstate2.acc.copy()
//     return (np.array([dda1,ddb1,ddg1,dda2,ddb2,ddg2]) @ np.array(dterms) + np.array([da1,db1,dg1,da2,db2,dg2]) @ np.array(ddterms) @ np.array([da1,db1,dg1,da2,db2,dg2]))

// def da1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def db1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*sinh(a2)*cos(b2) - cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def da2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*sinh(a2) - cosh(a1)*cos(b1)*cosh(a2)*cos(b2) - cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

// def db2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cos(b1)*sinh(a2)*sin(b2) - cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// def dg2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def db1D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*cos(b1)*sinh(a2)*cos(b2) + sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg1D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def db2D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// def dg2D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def dg1D12g1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg2D12g1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)
// # For the remaining nine functions of the upper triangular matrix use:
// # da2D12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)

// # da2D12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
// # db2D12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)

// # da2D12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
// # db2D12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
// # dg2D12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

// # db2D12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
// # dg2D12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

// # dg2D12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)

// # Third Derivatives of distance function (Only needed for implicit methods for jacobian construction)

// # da1D12a1 terms

// # def da1D12a1(state1, state2):
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def da1da1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*cosh(a2) - cosh(a1)*cos(b1)*sinh(a2)*cos(b2) - cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def db1da1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*sin(b1)*sinh(a2)*cos(b2) - sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg1da1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// # def da2da1D12a1(state1, state2): # same as da1da1D12a1 with arg swap
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return cosh(a1)*sinh(a2) - sinh(a1)*cos(b1)*cosh(a2)*cos(b2) - sinh(a1)*sin(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

// # def db2da1D12a1(state1, state2): # same as db1da1D12a1 with arg swap
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return sinh(a1)*cos(b1)*sinh(a2)*sin(b2) - sinh(a1)*sin(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// # def dg2da1D12a1(state1, state2): # same as dg1da1D12a1 with arg swap
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// # db1D12a1 terms

// # def db1D12a1(state1, state2):
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return cosh(a1)*sin(b1)*sinh(a2)*cos(b2) - cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def da1db1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*sin(b1)*sinh(a2)*cos(b2) - sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def db1db1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cos(b1)*sinh(a2)*cos(b2) + cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg1db1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def da2db1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*cosh(a2)*cos(b2) - cosh(a1)*cos(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

// def db2db1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return - cosh(a1)*sin(b1)*sinh(a2)*sin(b2) - cosh(a1)*cos(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// def dg2db1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return - cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// # dg1D12a1 terms

// # def dg1D12a1(state1, state2):
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def da1dg1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def db1dg1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def dg1dg1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def da2dg1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*sin(g1 - g2)

// def db2dg1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*sin(g1 - g2)

// def dg2dg1D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return - cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// # da2D12a1 terms

// # def da2D12a1(state1, state2):
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return sinh(a1)*sinh(a2) - cosh(a1)*cos(b1)*cosh(a2)*cos(b2) - cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

// def da1da2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sinh(a2) - sinh(a1)*cos(b1)*cosh(a2)*cos(b2) - sinh(a1)*sin(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

// def db1da2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*cosh(a2)*cos(b2) - cosh(a1)*cos(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

// def dg1da2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*sin(g1 - g2)

// def da2da2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*cosh(a2) - cosh(a1)*cos(b1)*sinh(a2)*cos(b2) - cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def db2da2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cos(b1)*cosh(a2)*sin(b2) - cosh(a1)*sin(b1)*cosh(a2)*cos(b2)*cos(g1 - g2)

// def dg2da2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return - cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*sin(g1 - g2)

// # db2D12a1 terms

// # def db2D12a1(state1, state2):
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return cosh(a1)*cos(b1)*sinh(a2)*sin(b2) - cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// def da1db2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*cos(b1)*sinh(a2)*sin(b2) - sinh(a1)*sin(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// def db1db2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return - cosh(a1)*sin(b1)*sinh(a2)*sin(b2) - cosh(a1)*cos(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// def dg1db2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*sin(g1 - g2)

// def da2db2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cos(b1)*cosh(a2)*sin(b2) - cosh(a1)*sin(b1)*cosh(a2)*cos(b2)*cos(g1 - g2)

// def db2db2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*cos(b1)*sinh(a2)*cos(b2) + cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg2db2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return - cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*sin(g1 - g2)

// # dg2D12a1 terms

// # def dg2D12a1(state1, state2):
// #     a1,b1,g1 = state1.pos.copy()
// #     a2,b2,g2 = state2.pos.copy()
// #     return -cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def da1dg2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def db1dg2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def dg1dg2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def da2dg2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*sin(g1 - g2)

// def db2dg2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*sin(g1 - g2)

// def dg2dg2D12a1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// # db1D12b1 terms



// def db1D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*cos(b1)*sinh(a2)*cos(b2) + sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg1D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def db2D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

// def dg2D12b1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

// def dg1D12g1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

// def dg2D12g1(state1, state2):
//     a1,b1,g1 = state1.pos.copy()
//     a2,b2,g2 = state2.pos.copy()
//     return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)
// # For the remaining nine functions of the upper triangular matrix use:
// # da2D12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)

// # da2D12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
// # db2D12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)

// # da2D12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
// # db2D12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
// # dg2D12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

// # db2D12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
// # dg2D12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

// # dg2D12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)



// # Derivative of metric components for pairwise coupling components
// # Derivative at state1 in terms of coordinates of both vertices
// # Needed for calculation of spring term contribution
// def dMetricTerms(state1):
//     a1,b1,g1 = state1.pos.copy()

// # Translational parameterization

//     # da1g11 = 0.
//     # db1g11 = 0.
//     # dg1g11 = 0.
    
//     # da2g11 = 0.
//     # db2g11 = 0.
//     # dg2g11 = 0.

//     # da1g22 = np.sinh(2. * a1)
//     # db1g22 = 0.
//     # dg1g22 = 0.

//     # da2g22 = 0
//     # db2g22 = 0.
//     # dg2g22 = 0.

//     # da1g33 = np.sinh(2. * a1) * np.cosh(b1)**2.
//     # db1g33 = np.sinh(2. * b1) * np.cosh(a1)**2.
//     # dg1g33 = 0.

//     # da2g33 = 0
//     # db2g33 = 0.
//     # dg2g33 = 0.

// # Rotational parameterization
    
//     da1g11 = 0.
//     db1g11 = 0.
//     dg1g11 = 0.
    
//     da2g11 = 0.
//     db2g11 = 0.
//     dg2g11 = 0.

//     da1g22 = np.sinh(2. * a1)
//     db1g22 = 0.
//     dg1g22 = 0.

//     da2g22 = 0
//     db2g22 = 0.
//     dg2g22 = 0.

//     da1g33 = np.sinh(2. * a1) * np.sin(b1)**2.
//     db1g33 = np.sin(2. * b1) * np.sinh(a1)**2.
//     dg1g33 = 0.

//     da2g33 = 0
//     db2g33 = 0.
//     dg2g33 = 0.

//     return np.array([
//         [da1g11,db1g11,dg1g11,da2g11,db2g11,dg2g11],
//         [da1g22,db1g22,dg1g22,da2g22,db2g22,dg2g22],
//         [da1g33,db1g33,dg1g33,da2g33,db2g33,dg2g33]
//     ])

// # Expression to generate coupling potential terms in system of odes and jacobian

// # First derivative term of coupling potential
// def da1V12(m, f, k, l, d12,  da1d12):
//     return k*(arccosh(d12) - l)*da1d12/(m * f * sqrt(d12**2. - 1.))

// # Second derivative term of coupling potential
// def da2da1V12(m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1):
//     # negative sign here so that the term can be added to geo terms later
//     return -k/(m*f*sqrt( d12**2. - 1. ))*( (da1d12*da2d12)/sqrt( d12**2. - 1.) + ( arccosh(d12) - l )*( da2d12da1 - da1d12*(da2f/f + d12*da2d12/(d12**2. - 1.)) ) )

// # Expression to generate rigidity constraint terms in system of odes and jacobian

// # Pairwise rigidity constraint
// def con12(l, d12):
//     return (arccosh(d12) - l)

// def dtcon12(l, d12,  dtd12):
//     return dtd12/(sqrt(d12**2. - 1.))

// def ddtcon12(l, d12,  dtd12, dttd12):
//     return 1./(sqrt(d12**2. - 1.)) * (dttd12 - d12*dtd12**2./(d12**2. - 1.))

// # Use first and second derivative functions from spring data above since constant becomes zero with derivative

// # First derivative term of rigidity constraint (for use in system of odes)
// def da1con12(m, f, lam, d12,  da1d12):
//     return lam*da1d12/(m * f * sqrt(d12**2. - 1.))

// # First derivative of velocity constriant for jacobian
// def da1dtcon12array(dstate1, dstate2, d12, dtd12, dterms, da1dterms):
//     # Replace dterms with da1dterms - i.e. derivative of dterms with respect to da1 typically given by the hessian
//     da1,db1,dg1 = dstate1.vel.copy()
//     da2,db2,dg2 = dstate2.vel.copy()
//     dterms = np.array(dterms)
//     da1dtd12arr = np.array([da1,db1,dg1,da2,db2,dg2]) @ np.array(da1dterms)
//     dda1dtd12arr = np.eye(6) @ np.array(dterms)
//     da1dtconterms  = (da1dtd12arr - (d12 * dterms * dtd12)/(d12**2. - 1))/np.sqrt(d12**2. - 1.)
//     dda1dtconterms = (dda1dtd12arr)/np.sqrt(d12**2. - 1.)
//     return np.array([
//         da1dtconterms[0], da1dtconterms[1], da1dtconterms[2],
//         dda1dtconterms[0],dda1dtconterms[1],dda1dtconterms[2],
//         da1dtconterms[3], da1dtconterms[4], da1dtconterms[5],
//         dda1dtconterms[3],dda1dtconterms[4],dda1dtconterms[5],
//     ])

// # First derivative of acceleration constraint for jacobian (WRONG need to consider derivative of dda1 terms...)
// # def da1ddtcon12(dstate1, dstate2, dterms, da1ddterms):
// #     da1,db1,dg1 = dstate1.vel.copy()
// #     da2,db2,dg2 = dstate2.vel.copy()
// #     dda1,ddb1,ddg1 = dstate1.acc.copy()
// #     dda2,ddb2,ddg2 = dstate2.acc.copy()
// #     return (np.array([dda1,ddb1,ddg1,dda2,ddb2,ddg2]) @ np.array(dterms) + np.array([da1,db1,dg1,da2,db2,dg2]) @ np.array(da1ddterms) @ np.array([da1,db1,dg1,da2,db2,dg2]))


// # Second derivative term of rigidity constraint (for use in jacobian)
// def da2da1con12(m, f, lam, d12, da1d12, da2d12, da2f, da2d12da1):
//     # negative sign here so that the term can be added to geo terms later
//     return -lam/(m*f*sqrt( d12**2. - 1. ))*( (da1d12*da2d12)/(d12**2. - 1.) + da1d12*da2f/f - da2d12da1 )

// -------------------------------------------------------------
// Geometry Information
// -------------------------------------------------------------

// Eigen::Matrix3d eucMetricTensor(Eigen::Vector3d){
//     Eigen::Matrix3d mat;
//     mat.setIdentity();
//     return mat;
// }

// Eigen::Vector3d eucChristoffel(State state){
//     Eigen::Vector3d vec;
//     vec.setZero();
//     return vec;
// }

// double eucDistance(Eigen::Vector3d pos1, Eigen::Vector3d pos2){
//     return std::sqrt((pos2 - pos1).dot(pos2 - pos1));
// }

// Geometry eucSpace = Geometry(
//     eucMetricTensor,
//     eucChristoffel,
//     eucDistance
//     );

// def eucMetricTensor(pos):
//     return np.identity(3)

// def eucChristoffel(state):
//     return np.zeros(3)

// def eucDistance(pos1, pos2):
//     return np.sqrt(np.dot(np.subtract(pos1.copy(),pos2.copy()),np.subtract(pos1.copy(),pos2.copy())))


// eucSpace = Geometry(
//     eucMetricTensor,
//     eucChristoffel,
//     eucDistance
//     )

// -------------------------------------------------------------
// Model Information
// -------------------------------------------------------------

// Eigen::Vector3d identityR3(Eigen::Vector3d coords){
//     return coords;
// }

// double unitScaling(Eigen::Vector3d pos){
//     return 1.;
// }

// Model eucModel = Model(identityR3,unitScaling);

// def identityR3(coords):
//     return coords

// def unitScaling(pos):
//     return 1.

// eucModel = Model(identityR3,unitScaling)

// -------------------------------------------------------------
// Obstacle/Bounding Ball Information
// -------------------------------------------------------------

// // Sphere Bounding Box

// //Default bounding box size (radius)
// double R = 6.;

// double distToSphere(Eigen::Vector3d pos){
//     return R - pos.norm();
// }

// Obstacle sphereObstacle = Obstacle(
//     distToSphere,
//     R
// );

// // def distToSphere(pos):
// //     return R-np.sqrt(pos[0]**2. + pos[1]**2. + pos[2]**2.)

// // sphereObstacle = Obstacle(
// //     distToSphere,
// //     R
// // )

// // Cube Bounding Box

// //Default bounding box size (cube side length)
// double boxSize = 6.;

// double distToBox(Eigen::Vector3d pos){
//     double xWall = boxSize - std::abs(pos[0]);
//     double yWall = boxSize - std::abs(pos[1]);
//     double zWall = boxSize - std::abs(pos[2]);

//     return std::min(xWall,std::min(yWall,zWall));
// }

// Obstacle boxObstacle = Obstacle(
//     distToBox,
//     boxSize
// );

// def distToBox(pos):
//     xWall = boxSize - abs(pos.x)
//     yWall = boxSize - abs(pos.y)
//     zWall = boxSize - abs(pos.z)

//     return min(xWall,min(yWall,zWall))

// boxObstacle = Obstacle(
//     distToBox,
//     boxSize
// )

//package stuff up for export
// AmbientSpace euclidean = AmbientSpace( eucSpace, eucModel, sphereObstacle);

#endif