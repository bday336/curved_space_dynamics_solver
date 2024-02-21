# // a class for integrating equations of motion called "integrator" and one specific implementation, Runge Kutta
# //derive is a function taking a state to state (now storing velocity and acceleration instead of position and velocity)
# //items fed into RungeKutta need to have the following methods available:
# //.add, .multiplyScalar, .clone

import numpy as np
from numpy import sin, cos, tan, cosh, sinh, tanh, arccosh, sqrt, array, zeros

from src.Computation.DataList import DataList
from src.Computation.State import State
from src.Computation.dState import dState

# //implementing the Rk4 Scheme for arbitrary classes that have clone add and multiplyScalar
# //will use this possibly on individual states, or on entire DataLists!
class Gauss1:
    """
    A class used to perform numerical integration via Explicit 4-stage Runge-Kutta method (RK4)

    ...

    Attributes
    ----------
    ambientSpace : object
        AmbientSpace object characterizing ambient space

    stepSize : float
        Integration step size

    Methods
    -------
    deriveFunc(dataList)
        Generates DataList consisting of dState characterizing dynamics acting on system (i.e. covariant acceleration (Christoffel symbol terms), external/internal potentials, etc.)
        Returns DataList of dStates

    step(dataList)
        Perform one integration step via the RK4 algorithm on system described by dataList
        Returns updated clone of dataList
    """

    def __init__(self, ambientSpace, stepSize):
        # self.derive=derive
        self.ambientSpace = ambientSpace
        self.stepSize = stepSize

        self.ks = 1.
        self.x = 1.
        self.m = 1.

    def arrayToDataList(self, arr, dataList):
        # Array is an np.array
        temparr = []

        # Contributions from ambient space geometry
        for a in range(len(dataList.data)):
            # print(a)
            temparr.append(State(arr[a*2*3:a*2*3 + 3],arr[a*2*3 + 3:a*2*3 + 6]))

        return DataList(temparr, connectivity=dataList.connectivity)

    def dynfunc(self, dataList, params = None):
        temparr = []

        # Contributions from ambient space geometry
        for a in range(len(dataList.data)):
            # print(a)
            temparr.append(self.ambientSpace.acceleration(dataList.data[a].clone()))

        # ks,x = [100,1]
        # m = 1.

        # Contributions from coupling potential
        if len(dataList.connectivity) != 0:
            for b in dataList.connectivity:
                # print(b)
                b.sort()
                # print(b)
                # Distance Function
                d12 = self.ambientSpace.geometry.funcDict["d12"](dataList.data[b[0]],dataList.data[b[1]])
                # First derivatives of distance function
                da1d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                db1d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[0]],dataList.data[b[1]])
                dg1d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[0]],dataList.data[b[1]])
                da2d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[1]],dataList.data[b[0]])
                db2d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[1]],dataList.data[b[0]])
                dg2d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[1]],dataList.data[b[0]])

                a1,b1,g1 = dataList.data[b[0]].pos.copy()
                a2,b2,g2 = dataList.data[b[1]].pos.copy()

                spa1 = (self.ks*(arccosh(d12) - self.x)*da1d12)/(self.m*sqrt(d12**2. - 1.))
                spb1 = (self.ks*(arccosh(d12) - self.x)*db1d12)/(self.m*sinh(a1)**2. * sqrt(d12**2. - 1.))
                spg1 = (self.ks*(arccosh(d12) - self.x)*dg1d12)/(self.m*sinh(a1)**2. * sin(b1)**2. * sqrt(d12**2. - 1.))
                spa2 = (self.ks*(arccosh(d12) - self.x)*da2d12)/(self.m*sqrt(d12**2. - 1.))
                spb2 = (self.ks*(arccosh(d12) - self.x)*db2d12)/(self.m*sinh(a2)**2. * sqrt(d12**2. - 1.))
                spg2 = (self.ks*(arccosh(d12) - self.x)*dg2d12)/(self.m*sinh(a2)**2. * sin(b2)**2. * sqrt(d12**2. - 1.))

                spdStates = [dState(zeros(3),array([spa1,spb1,spg1])), dState(zeros(3),array([spa2,spb2,spg2]))]

                temparr[b[0]] = temparr[b[0]].sub(spdStates[0])
                temparr[b[1]] = temparr[b[1]].sub(spdStates[1])

        # res_array = []
        # for c in temparr:
        #     res_array.append(c.vel.copy())
        #     res_array.append(c.acc.copy())

        # res_array = np.array(res_array).flatten

        # return [DataList(temparr, connectivity=dataList.connectivity), res_array]
        return DataList(temparr, connectivity=dataList.connectivity)

    # Jacobian will be constructed to be largely diagonal to help with efficiently for larger (sparse matrix) systems - block diagonal jacobian
    # This is in contrast with previous iteration of the solver design where top half was velocity residual expressions and bottom half were acceleration residual expressions
    def dynjac(self, dataList, params = None):
        # Dimension of ambient space (hard coded to 3D at the moment)
        ambient_dim = 3
        # Number of vertices in system
        vert_num = len(dataList.data)
        # Initialize jacobian matrix (multiply by 2 since considering residuals of velocity and accelerations for each vertex)
        jacobian = np.zeros((2*ambient_dim*vert_num,2*ambient_dim*vert_num))

        # Populate vertex data (block diagonal information in jacobian)
        for a in range(vert_num):
            tempmat = self.ambientSpace.geometry.funcDict["vertex_jac"](dataList.data[a])
            # print(a*2*ambient_dim )
            # print(a*2*ambient_dim + 2 * ambient_dim)
            jacobian[a*2*ambient_dim : a*2*ambient_dim + 2 * ambient_dim , a*2*ambient_dim : a*2*ambient_dim + 2 * ambient_dim] = tempmat.copy()

        # ks,x = [1,1]
        # m = 1.

        # Contributions from coupling potential (right now only pairwise coupling)
        if len(dataList.connectivity) != 0:
            for b in dataList.connectivity:
                b.sort()
                # Distance Function
                d12 = self.ambientSpace.geometry.funcDict["d12"](dataList.data[b[0]],dataList.data[b[1]])
                # First derivatives of distance function
                da1d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                db1d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[0]],dataList.data[b[1]])
                dg1d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[0]],dataList.data[b[1]])
                da2d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[1]],dataList.data[b[0]])
                db2d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[1]],dataList.data[b[0]])
                dg2d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[1]],dataList.data[b[0]])
                # Second derivatives of distance function
                da1d12a1 = self.ambientSpace.geometry.funcDict["da1da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                db1d12a1 = self.ambientSpace.geometry.funcDict["db1da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                dg1d12a1 = self.ambientSpace.geometry.funcDict["dg1da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                da2d12a1 = self.ambientSpace.geometry.funcDict["da2da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                db2d12a1 = self.ambientSpace.geometry.funcDict["db2da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                dg2d12a1 = self.ambientSpace.geometry.funcDict["dg2da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                
                da1d12b1 = db1d12a1
                db1d12b1 = self.ambientSpace.geometry.funcDict["db1db1d12"](dataList.data[b[0]],dataList.data[b[1]])
                dg1d12b1 = self.ambientSpace.geometry.funcDict["dg1db1d12"](dataList.data[b[0]],dataList.data[b[1]])
                da2d12b1 = self.ambientSpace.geometry.funcDict["db2da1d12"](dataList.data[b[1]],dataList.data[b[0]])
                db2d12b1 = self.ambientSpace.geometry.funcDict["db2db1d12"](dataList.data[b[0]],dataList.data[b[1]])
                dg2d12b1 = self.ambientSpace.geometry.funcDict["dg2db1d12"](dataList.data[b[0]],dataList.data[b[1]])

                da1d12g1 = dg1d12a1
                db1d12g1 = dg1d12b1
                dg1d12g1 = self.ambientSpace.geometry.funcDict["dg1dg1d12"](dataList.data[b[0]],dataList.data[b[1]])
                da2d12g1 = self.ambientSpace.geometry.funcDict["dg2da1d12"](dataList.data[b[1]],dataList.data[b[0]])
                db2d12g1 = self.ambientSpace.geometry.funcDict["dg2db1d12"](dataList.data[b[1]],dataList.data[b[0]])
                dg2d12g1 = self.ambientSpace.geometry.funcDict["dg2dg1d12"](dataList.data[b[0]],dataList.data[b[1]])

                da1d12a2 = da2d12a1
                db1d12a2 = da2d12b1
                dg1d12a2 = da2d12g1
                da2d12a2 = self.ambientSpace.geometry.funcDict["da1da1d12"](dataList.data[b[1]],dataList.data[b[0]])
                db2d12a2 = self.ambientSpace.geometry.funcDict["db1da1d12"](dataList.data[b[1]],dataList.data[b[0]])
                dg2d12a2 = self.ambientSpace.geometry.funcDict["dg1da1d12"](dataList.data[b[1]],dataList.data[b[0]])

                da1d12b2 = db2d12a1
                db1d12b2 = db2d12b1
                dg1d12b2 = db2d12g1
                da2d12b2 = db2d12a2
                db2d12b2 = self.ambientSpace.geometry.funcDict["db1db1d12"](dataList.data[b[1]],dataList.data[b[0]])
                dg2d12b2 = self.ambientSpace.geometry.funcDict["dg1db1d12"](dataList.data[b[1]],dataList.data[b[0]])

                da1d12g2 = dg2d12a1
                db1d12g2 = dg2d12b1
                dg1d12g2 = dg2d12g1
                da2d12g2 = dg2d12a2
                db2d12g2 = dg2d12b2
                dg2d12g2 = self.ambientSpace.geometry.funcDict["dg1dg1d12"](dataList.data[b[1]],dataList.data[b[0]])

                a1,b1,g1 = dataList.data[b[0]].pos.copy()
                a2,b2,g2 = dataList.data[b[1]].pos.copy()

                #---------- particle 1

                da1da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da1d12, da1d12, 0., da1d12a1)
                db1da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da1d12, db1d12, 0., db1d12a1)
                dg1da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da1d12, dg1d12, 0., dg1d12a1)

                da2da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da1d12, da2d12, 0., da2d12a1)
                db2da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da1d12, db2d12, 0., db2d12a1)
                dg2da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da1d12, dg2d12, 0., dg2d12a1)

                #----------

                da1db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2., self.ks, self.x, d12, db1d12, da1d12, sinh(2.*a1), da1d12b1)
                db1db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2., self.ks, self.x, d12, db1d12, db1d12, 0.,          db1d12b1)
                dg1db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2., self.ks, self.x, d12, db1d12, dg1d12, 0.,          dg1d12b1)

                da2db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2., self.ks, self.x, d12, db1d12, da2d12, 0.,          da2d12b1)
                db2db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2., self.ks, self.x, d12, db1d12, db2d12, 0.,          db2d12b1)
                dg2db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2., self.ks, self.x, d12, db1d12, dg2d12, 0.,          dg2d12b1)

                #------------

                da1dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2.*sin(b1)**2., self.ks, self.x, d12, dg1d12, da1d12, sinh(2.*a1)*sin(b1)**2., da1d12g1)
                db1dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2.*sin(b1)**2., self.ks, self.x, d12, dg1d12, db1d12, sin(2.*b1)*sinh(a1)**2., db1d12g1)
                dg1dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2.*sin(b1)**2., self.ks, self.x, d12, dg1d12, dg1d12, 0.,                      dg1d12g1)

                da2dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2.*sin(b1)**2., self.ks, self.x, d12, dg1d12, da2d12, 0.,                      da2d12g1)
                db2dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2.*sin(b1)**2., self.ks, self.x, d12, dg1d12, db2d12, 0.,                      db2d12g1)
                dg2dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a1)**2.*sin(b1)**2., self.ks, self.x, d12, dg1d12, dg2d12, 0.,                      dg2d12g1)

                # Contribution to particle 1 dynamics (block diagonal contribution)
                part1_geoterms = np.array([
                    [da1da1V12,db1da1V12,dg1da1V12,  0., 0., 0.],
                    [da1db1V12,db1db1V12,dg1db1V12,  0., 0., 0.],
                    [da1dg1V12,db1dg1V12,dg1dg1V12,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.]
                ])

                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] + part1_geoterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part1_spterms = np.array([
                    [da2da1V12,db2da1V12,dg2da1V12, 0., 0., 0.],
                    [da2db1V12,db2db1V12,dg2db1V12, 0., 0., 0.],
                    [da2dg1V12,db2dg1V12,dg2dg1V12, 0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.]
                ])

                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part1_spterms.copy()

                #------------- particle 2

                da1da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da2d12, da1d12, 0., da1d12a2)
                db1da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da2d12, db1d12, 0., db1d12a2)
                dg1da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da2d12, dg1d12, 0., dg1d12a2)

                da2da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da2d12, da2d12, 0., da2d12a2)
                db2da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da2d12, db2d12, 0., db2d12a2)
                dg2da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, 1., self.ks, self.x, d12, da2d12, dg2d12, 0., dg2d12a2)

                #----------

                da1db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2., self.ks, self.x, d12, db2d12, da1d12, 0.,          da1d12b2)
                db1db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2., self.ks, self.x, d12, db2d12, db1d12, 0.,          db1d12b2)
                dg1db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2., self.ks, self.x, d12, db2d12, dg1d12, 0.,          dg1d12b2)

                da2db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2., self.ks, self.x, d12, db2d12, da2d12, sinh(2.*a2), da2d12b2)
                db2db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2., self.ks, self.x, d12, db2d12, db2d12, 0.,          db2d12b2)
                dg2db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2., self.ks, self.x, d12, db2d12, dg2d12, 0.,          dg2d12b2)

                #------------

                da1dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2.*sin(b2)**2., self.ks, self.x, d12, dg2d12, da1d12, 0.,                      da1d12g2)
                db1dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2.*sin(b2)**2., self.ks, self.x, d12, dg2d12, db1d12, 0.,                      db1d12g2)
                dg1dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2.*sin(b2)**2., self.ks, self.x, d12, dg2d12, dg1d12, 0.,                      dg1d12g2)

                da2dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2.*sin(b2)**2., self.ks, self.x, d12, dg2d12, da2d12, sinh(2.*a2)*sin(b2)**2., da2d12g2)
                db2dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2.*sin(b2)**2., self.ks, self.x, d12, dg2d12, db2d12, sin(2.*b2)*sinh(a2)**2., db2d12g2)
                dg2dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative"](self.m, sinh(a2)**2.*sin(b2)**2., self.ks, self.x, d12, dg2d12, dg2d12, 0.,                      dg2d12g2)

                # Contribution to particle 1 dynamics (block diagonal contribution)
                part2_geoterms = np.array([
                    [da2da2V12,db2da2V12,dg2da2V12,  0., 0., 0.],
                    [da2db2V12,db2db2V12,dg2db2V12,  0., 0., 0.],
                    [da2dg2V12,db2dg2V12,dg2dg2V12,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.]
                ])

                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part2_geoterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part2_spterms = np.array([
                    [da1da2V12,db1da2V12,dg1da2V12, 0., 0., 0.],
                    [da1db2V12,db1db2V12,dg1db2V12, 0., 0., 0.],
                    [da1dg2V12,db1dg2V12,dg1dg2V12, 0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.]
                ])

                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] + part2_spterms.copy()

        # fullj = np.eye(2*ambient_dim*vert_num) - jacobian
        return jacobian

    def difffuncgauss1s(self, dataList, k1):
        # Set to run Gauss 1-stage method
        a11 = 1./2.

        return np.array([
            k1.toArray() - self.dynfunc(self.arrayToDataList(dataList.toArray() + (a11*k1.toArray())*self.stepSize, dataList)).toArray()
        ]).flatten()

    def step(self, dataList, tol = 1e-15, imax = 100):
        # Set to run Gauss 1-stage method
        a11 = 1./2.
        bs1 = 1.

        # Initial Guess - Explicit Euler
        k = self.dynfunc(dataList)
        x1guess = dataList.toArray() + (1./2.)*self.stepSize*k.toArray()
        k1 = self.dynfunc(self.arrayToDataList(x1guess, dataList))

        # Check Error Before iterations
        er = self.difffuncgauss1s(dataList, k1)

        # Begin Iterations
        for a in range(imax):
            if np.linalg.norm(er) >= tol:
                j1 = self.dynjac(self.arrayToDataList(dataList.toArray() + (a11*k1.toArray())*self.stepSize, dataList))
                
                fullj = np.block([
                    [np.eye(j1.shape[0]) - self.stepSize*a11*j1]
                ])

                linsolve = np.linalg.solve(fullj,-er)

                k1 = self.arrayToDataList(k1.toArray() + linsolve[0:j1.shape[0]],dataList)

                er = self.difffuncgauss1s(dataList, k1)
            else:
                break


        newarr = dataList.toArray() + self.stepSize*(bs1*k1.toArray())
        return self.arrayToDataList(newarr, dataList)
    
        # def difffuncgauss1s(self, dataList, k1DataList, params = None):
        #     # Set to run Gauss 1-stage method
        #     a11 = 1./2.

        #     return (k1DataList.sub(self.dynfunc(dataList.add(k1DataList.multiplyScalar(a11*self.stepSize)), params = None)))

    # def step(self, dataList, params = None, tol = 1e-15, imax = 100):
    #     # Set to run Gauss 1-stage method
    #     a11 = 1./2.
    #     bs1 = 1.

    #     # Initial Guess - Explicit Euler
    #     kDataList = self.dynfunc(dataList, params = None)
    #     x1guess = dataList.add(kDataList.multiplyScalar(.5 * self.stepSize))
    #     k1DataList = self.dynfunc(x1guess, params = None)

    #     # Check Error Before iterations
    #     # er = self.difffuncgauss1s(dataList, k1DataList, params = None)
    #     erDataList = k1DataList.sub(self.dynfunc(dataList.add(k1DataList.multiplyScalar(a11*self.stepSize)), params = None))
    #     res_array = []
    #     for c in erDataList.data:
    #         res_array.append(c.vel.copy())
    #         res_array.append(c.acc.copy())

    #     res_array = np.array(res_array).flatten

    #     # Begin Iterations
    #     for a in range(imax): #(er.norm() >= tol and counter <= imax):
    #         if erDataList.norm() >= tol:
    #             break
    #         else:
    #             j1 = self.dynjac(dataList.add(k1DataList.multiplyScalar(a11*self.stepSize)), params = None)
                
    #             fullj = np.block([
    #                 [np.eye(j1.shape[0]) - self.stepSize*a11*j1]
    #             ])

    #             linsolve = np.linalg.solve(fullj,-res_array)

    #             k1 = k1 + linsolve[0:j1.shape[0]]

    #             er = self.difffuncgauss1s(startvec, params, dynfunc, k1, dt)

    #     startvec = startvec + dt*(bs1*k1)
    #     return startvec.copy()

    # //step forwards one timestep
    # def step(self, dataList):

    #     # k1,k2,k3,k4;
    #     # temp;

    #     # //get the derivative
    #     k1 = self.deriveFunc(dataList)
    #     # print(k1.data[1].vel)
    #     temp1 = dataList.clone().add(k1.clone().multiplyScalar(self.stepSize/2.))
    #     # k1.multiplyScalar(self.ep)

    #     # //get k2
    #     # temp = dataList.clone().add(k1.clone().multiplyScalar(0.5))
    #     k2 = self.deriveFunc(temp1)
    #     temp2 = dataList.clone().add(k2.clone().multiplyScalar(self.stepSize/2.))
    #     # k2.multiplyScalar(self.stepSize)

    #     # //get k3
    #     # temp=dataList.clone().add(k2.clone().multiplyScalar(0.5))
    #     k3=self.deriveFunc(temp2)
    #     temp3 = dataList.clone().add(k2.clone().multiplyScalar(self.stepSize))
    #     # k3.multiplyScalar(self.stepSize)

    #     # //get k4
    #     # temp=dataList.clone().add(k3.multiplyScalar(1.))
    #     k4=self.deriveFunc(temp3)
    #     # k4.multiplyScalar(self.ep)

    #     # //add up results:
    #     total = k1 #//scale factor 1
    #     total.add(k2.multiplyScalar(2))
    #     total.add(k3.multiplyScalar(2))
    #     total.add(k4) #//scale factor 1
    #     total.multiplyScalar(self.stepSize/6.)

    #     # print(dataList.clone().updateBy(total).data[0].pos)
    #     # print(total.data[1].vel)
    #     # //move ahead one step
    #     return dataList.clone().updateBy(total)