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
class Gauss2test:
    """
    A class used to perform numerical integration via Implicit 2-stage Gauss method (GS2)

    ...

    Attributes
    ----------
    ambientSpace : object
        AmbientSpace object characterizing ambient space

    stepSize : float
        Integration step size

    Methods
    -------
    arrayToDataList(arr, dataList)
        Generates DataList consisting of list of States in arr with the same connectivity as dataList
        Returns DataList of States

    dynfunc(dataList)
        Generates DataList consisting of list of dStates corresponding with system of linear odes for simulation system
        Returns DataList of dStates

    dynjac(dataList)
        Generates jacobian matrix for solving the system of odes for simulation system with 2-step Gauss collocation
        Returns jacobian matrix (np.array)

    difffunc(dataList, k1, k2)
        Generates np.array of residuals of odes for simulation system with 1-step Gauss collocation using DataList objects dataList (initial condition for integration step), k1 (data at stage 1 of Gauss collocation), and k2 (data at stage 2 of Gauss collocation)
        Returns np.array

    step(dataList)
        Perform one integration step via the GS2 algorithm on system described by dataList
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

                spa1 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[0,0], self.ks, self.x, d12,  da1d12)
                spb1 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[1,1], self.ks, self.x, d12,  db1d12)
                spg1 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[2,2], self.ks, self.x, d12,  dg1d12)
                spa2 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[0,0], self.ks, self.x, d12,  da2d12)
                spb2 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[1,1], self.ks, self.x, d12,  db2d12)
                spg2 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[2,2], self.ks, self.x, d12,  dg2d12)


                # spa1 = (self.ks*(arccosh(d12) - self.x)*da1d12)/(self.m*sqrt(d12**2. - 1.))
                # spb1 = (self.ks*(arccosh(d12) - self.x)*db1d12)/(self.m*sinh(a1)**2. * sqrt(d12**2. - 1.))
                # spg1 = (self.ks*(arccosh(d12) - self.x)*dg1d12)/(self.m*sinh(a1)**2. * sin(b1)**2. * sqrt(d12**2. - 1.))
                # spa2 = (self.ks*(arccosh(d12) - self.x)*da2d12)/(self.m*sqrt(d12**2. - 1.))
                # spb2 = (self.ks*(arccosh(d12) - self.x)*db2d12)/(self.m*sinh(a2)**2. * sqrt(d12**2. - 1.))
                # spg2 = (self.ks*(arccosh(d12) - self.x)*dg2d12)/(self.m*sinh(a2)**2. * sin(b2)**2. * sqrt(d12**2. - 1.))

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

                metric1diag = self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos.copy()).diagonal()
                metric2diag = self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos.copy()).diagonal()

                dmetric_terms1 = self.ambientSpace.geometry.funcDict["dmetric_terms"](dataList.data[b[0]])
                dmetric_terms2 = self.ambientSpace.geometry.funcDict["dmetric_terms"](dataList.data[b[1]])

                # a1,b1,g1 = dataList.data[b[0]].pos.copy()
                # a2,b2,g2 = dataList.data[b[1]].pos.copy()

                #---------- particle 1

                da1da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[0], self.ks, self.x, d12, da1d12, da1d12, dmetric_terms1[0,0], da1d12a1)
                db1da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[0], self.ks, self.x, d12, da1d12, db1d12, dmetric_terms1[0,1], db1d12a1)
                dg1da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[0], self.ks, self.x, d12, da1d12, dg1d12, dmetric_terms1[0,2], dg1d12a1)

                da2da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[0], self.ks, self.x, d12, da1d12, da2d12, dmetric_terms1[0,3], da2d12a1)
                db2da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[0], self.ks, self.x, d12, da1d12, db2d12, dmetric_terms1[0,4], db2d12a1)
                dg2da1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[0], self.ks, self.x, d12, da1d12, dg2d12, dmetric_terms1[0,5], dg2d12a1)

                #----------

                da1db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[1], self.ks, self.x, d12, db1d12, da1d12, dmetric_terms1[1,0], da1d12b1)
                db1db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[1], self.ks, self.x, d12, db1d12, db1d12, dmetric_terms1[1,1], db1d12b1)
                dg1db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[1], self.ks, self.x, d12, db1d12, dg1d12, dmetric_terms1[1,2], dg1d12b1)

                da2db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[1], self.ks, self.x, d12, db1d12, da2d12, dmetric_terms1[1,3], da2d12b1)
                db2db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[1], self.ks, self.x, d12, db1d12, db2d12, dmetric_terms1[1,4], db2d12b1)
                dg2db1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[1], self.ks, self.x, d12, db1d12, dg2d12, dmetric_terms1[1,5], dg2d12b1)

                #------------

                da1dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[2], self.ks, self.x, d12, dg1d12, da1d12, dmetric_terms1[2,0], da1d12g1)
                db1dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[2], self.ks, self.x, d12, dg1d12, db1d12, dmetric_terms1[2,1], db1d12g1)
                dg1dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[2], self.ks, self.x, d12, dg1d12, dg1d12, dmetric_terms1[2,2], dg1d12g1)

                da2dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[2], self.ks, self.x, d12, dg1d12, da2d12, dmetric_terms1[2,3], da2d12g1)
                db2dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[2], self.ks, self.x, d12, dg1d12, db2d12, dmetric_terms1[2,4], db2d12g1)
                dg2dg1V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric1diag[2], self.ks, self.x, d12, dg1d12, dg2d12, dmetric_terms1[2,5], dg2d12g1)

                # Contribution to particle 1 dynamics (block diagonal contribution)
                part1_geoterms = np.array([
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [da1da1V12,db1da1V12,dg1da1V12,  0., 0., 0.],
                    [da1db1V12,db1db1V12,dg1db1V12,  0., 0., 0.],
                    [da1dg1V12,db1dg1V12,dg1dg1V12,  0., 0., 0.]
                ])

                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] + part1_geoterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part1_spterms = np.array([
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [da2da1V12,db2da1V12,dg2da1V12, 0., 0., 0.],
                    [da2db1V12,db2db1V12,dg2db1V12, 0., 0., 0.],
                    [da2dg1V12,db2dg1V12,dg2dg1V12, 0., 0., 0.]
                ])

                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part1_spterms.copy()

                #------------- particle 2

                da1da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[0], self.ks, self.x, d12, da2d12, da1d12, dmetric_terms2[0,0], da1d12a2)
                db1da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[0], self.ks, self.x, d12, da2d12, db1d12, dmetric_terms2[0,1], db1d12a2)
                dg1da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[0], self.ks, self.x, d12, da2d12, dg1d12, dmetric_terms2[0,2], dg1d12a2)

                da2da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[0], self.ks, self.x, d12, da2d12, da2d12, dmetric_terms2[0,3], da2d12a2)
                db2da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[0], self.ks, self.x, d12, da2d12, db2d12, dmetric_terms2[0,4], db2d12a2)
                dg2da2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[0], self.ks, self.x, d12, da2d12, dg2d12, dmetric_terms2[0,5], dg2d12a2)

                #----------

                da1db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[1], self.ks, self.x, d12, db2d12, da1d12, dmetric_terms2[1,0], da1d12b2)
                db1db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[1], self.ks, self.x, d12, db2d12, db1d12, dmetric_terms2[1,1], db1d12b2)
                dg1db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[1], self.ks, self.x, d12, db2d12, dg1d12, dmetric_terms2[1,2], dg1d12b2)

                da2db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[1], self.ks, self.x, d12, db2d12, da2d12, dmetric_terms2[1,3], da2d12b2)
                db2db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[1], self.ks, self.x, d12, db2d12, db2d12, dmetric_terms2[1,4], db2d12b2)
                dg2db2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[1], self.ks, self.x, d12, db2d12, dg2d12, dmetric_terms2[1,5], dg2d12b2)

                #------------

                da1dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[2], self.ks, self.x, d12, dg2d12, da1d12, dmetric_terms2[2,0], da1d12g2)
                db1dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[2], self.ks, self.x, d12, dg2d12, db1d12, dmetric_terms2[2,1], db1d12g2)
                dg1dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[2], self.ks, self.x, d12, dg2d12, dg1d12, dmetric_terms2[2,2], dg1d12g2)

                da2dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[2], self.ks, self.x, d12, dg2d12, da2d12, dmetric_terms2[2,3], da2d12g2)
                db2dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[2], self.ks, self.x, d12, dg2d12, db2d12, dmetric_terms2[2,4], db2d12g2)
                dg2dg2V12 = self.ambientSpace.geometry.funcDict["coupling_derivative2"](self.m, metric2diag[2], self.ks, self.x, d12, dg2d12, dg2d12, dmetric_terms2[2,5], dg2d12g2)

                # Contribution to particle 1 dynamics (block diagonal contribution)
                part2_geoterms = np.array([
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [da2da2V12,db2da2V12,dg2da2V12,  0., 0., 0.],
                    [da2db2V12,db2db2V12,dg2db2V12,  0., 0., 0.],
                    [da2dg2V12,db2dg2V12,dg2dg2V12,  0., 0., 0.]
                ])

                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part2_geoterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part2_spterms = np.array([
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [0., 0., 0.,  0., 0., 0.],
                    [da1da2V12,db1da2V12,dg1da2V12, 0., 0., 0.],
                    [da1db2V12,db1db2V12,dg1db2V12, 0., 0., 0.],
                    [da1dg2V12,db1dg2V12,dg1dg2V12, 0., 0., 0.]
                ])

                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] + part2_spterms.copy()

        # fullj = np.eye(2*ambient_dim*vert_num) - jacobian
        return jacobian

    def difffunc(self, dataList, k1, k2):
        # Set to run Gauss 2-stage method
        a11,a12 = [1./4., 1./4. - np.sqrt(3.)/6.]
        a21,a22 = [1./4. + np.sqrt(3.)/6., 1./4.]

        return np.array([
            k1.toArray() - self.dynfunc(self.arrayToDataList(dataList.toArray() + (a11*k1.toArray() + a12*k2.toArray())*self.stepSize, dataList)).toArray(),
            k2.toArray() - self.dynfunc(self.arrayToDataList(dataList.toArray() + (a21*k1.toArray() + a22*k2.toArray())*self.stepSize, dataList)).toArray()
        ]).flatten()

    def step(self, dataList, tol = 1e-15, imax = 100):
        # Set to run Gauss 2-stage method
        a11,a12 = [1./4., 1./4. - np.sqrt(3.)/6.]
        a21,a22 = [1./4. + np.sqrt(3.)/6., 1./4.]
        bs1,bs2 = [1./2., 1./2.]

        # Initial Guess - Explicit Euler
        k = self.dynfunc(dataList)
        x1guess = dataList.toArray() + (1./2. - np.sqrt(3.)/6.)*self.stepSize*k.toArray()
        x2guess = dataList.toArray() + (1./2. + np.sqrt(3.)/6.)*self.stepSize*k.toArray()
        k1 = self.dynfunc(self.arrayToDataList(x1guess, dataList))
        k2 = self.dynfunc(self.arrayToDataList(x2guess, dataList))

        # Check Error Before iterations
        er = self.difffunc(dataList, k1, k2)

        # Begin Iterations
        for a in range(imax):
            if np.linalg.norm(er) >= tol:
                j1 = self.dynjac(self.arrayToDataList(dataList.toArray() + (a11*k1.toArray() + a12*k2.toArray())*self.stepSize, dataList))
                j2 = self.dynjac(self.arrayToDataList(dataList.toArray() + (a21*k1.toArray() + a22*k2.toArray())*self.stepSize, dataList))
                
                fullj = np.block([
                    [np.eye(j1.shape[0]) - self.stepSize*a11*j1, -self.stepSize*a12*j1],
                    [-self.stepSize*a21*j2, np.eye(j2.shape[0]) - self.stepSize*a22*j2]
                ])

                linsolve = np.linalg.solve(fullj,-er)

                k1 = self.arrayToDataList(k1.toArray() + linsolve[0:j1.shape[0]],dataList)
                k2 = self.arrayToDataList(k2.toArray() + linsolve[j2.shape[0]:2*j2.shape[0]],dataList)

                er = self.difffunc(dataList, k1, k2)
            else:
                break

        np.set_printoptions(precision=1)
        print(fullj)


        newarr = dataList.toArray() + self.stepSize*(bs1*k1.toArray() + bs2*k2.toArray())
        return self.arrayToDataList(newarr, dataList)