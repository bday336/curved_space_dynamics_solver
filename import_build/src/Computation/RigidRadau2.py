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
class RigidRadau2:
    """
    A class used to perform numerical integration via Implicit 2-stage Radau method with rigidity constraints(RigidRS2)

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
        Generates jacobian matrix for solving the system of odes for simulation system with 2-step Radau collocation
        Returns jacobian matrix (np.array)

    difffunc(dataList, k1, k2)
        Generates np.array of residuals of odes for simulation system with 1-step Radau collocation using DataList objects dataList (initial condition for integration step), k1 (data at stage 1 of Radau collocation), and k2 (data at stage 2 of Radau collocation)
        Returns np.array

    step(dataList)
        Perform one integration step via the RS2 algorithm on system described by dataList
        Returns updated clone of dataList
    """

    def __init__(self, ambientSpace, stepSize):
        # self.derive=derive
        self.ambientSpace = ambientSpace
        self.stepSize = stepSize

        # self.ks = 0.
        self.x = 1.
        self.m = 1.

    def arrayToDataList(self, arr, dataList):
        # Array is an np.array
        temparr = []

        for a in range(len(dataList.data)):
            # print(a)
            temparr.append(State(arr[a*2*3:a*2*3 + 3],arr[a*2*3 + 3:a*2*3 + 6]))

        return DataList(temparr, connectivity=dataList.connectivity, rig_connectivity=dataList.rig_connectivity)

    def dynfunc(self, dataList, stage, params = None):
        temparr = []
        tempconarr = []
        # print(stage)

        # Contributions from ambient space geometry
        for a in range(len(dataList.data)):
            # print(a)
            temparr.append(self.ambientSpace.acceleration(dataList.data[a].clone()))

        # ks,x = [100,1]
        # m = 1.

        # # Contributions from coupling potential
        # if len(dataList.connectivity) != 0:
        #     for b in dataList.connectivity:
        #         # print(b)
        #         b.sort()
        #         # print(b)
        #         # Distance Function
        #         d12 = self.ambientSpace.geometry.funcDict["d12"](dataList.data[b[0]],dataList.data[b[1]])
        #         # First derivatives of distance function
        #         da1d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[0]],dataList.data[b[1]])
        #         db1d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[0]],dataList.data[b[1]])
        #         dg1d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[0]],dataList.data[b[1]])
        #         da2d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[1]],dataList.data[b[0]])
        #         db2d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[1]],dataList.data[b[0]])
        #         dg2d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[1]],dataList.data[b[0]])

        #         a1,b1,g1 = dataList.data[b[0]].pos.copy()
        #         a2,b2,g2 = dataList.data[b[1]].pos.copy()

        #         spa1 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[0,0], self.ks, self.x, d12,  da1d12)
        #         spb1 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[1,1], self.ks, self.x, d12,  db1d12)
        #         spg1 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[2,2], self.ks, self.x, d12,  dg1d12)
        #         spa2 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[0,0], self.ks, self.x, d12,  da2d12)
        #         spb2 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[1,1], self.ks, self.x, d12,  db2d12)
        #         spg2 = self.ambientSpace.geometry.funcDict["coupling_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[2,2], self.ks, self.x, d12,  dg2d12)


        #         # spa1 = (self.ks*(arccosh(d12) - self.x)*da1d12)/(self.m*sqrt(d12**2. - 1.))
        #         # spb1 = (self.ks*(arccosh(d12) - self.x)*db1d12)/(self.m*sinh(a1)**2. * sqrt(d12**2. - 1.))
        #         # spg1 = (self.ks*(arccosh(d12) - self.x)*dg1d12)/(self.m*sinh(a1)**2. * sin(b1)**2. * sqrt(d12**2. - 1.))
        #         # spa2 = (self.ks*(arccosh(d12) - self.x)*da2d12)/(self.m*sqrt(d12**2. - 1.))
        #         # spb2 = (self.ks*(arccosh(d12) - self.x)*db2d12)/(self.m*sinh(a2)**2. * sqrt(d12**2. - 1.))
        #         # spg2 = (self.ks*(arccosh(d12) - self.x)*dg2d12)/(self.m*sinh(a2)**2. * sin(b2)**2. * sqrt(d12**2. - 1.))

        #         spdStates = [dState(zeros(3),array([spa1,spb1,spg1])), dState(zeros(3),array([spa2,spb2,spg2]))]

        # Contributions from rigidity constriants
        if len(dataList.rig_connectivity) != 0:
            for b in dataList.rig_connectivity:
                # print(b)
                # Change format to be in the form of 2-stage radau method
                if type(b[-1]).__name__ != 'list':
                    # print("reformat")
                    b[-1] = [b[-1],b[-1]]
                # print(b)
                indx = b[0:2]
                indx.sort()
                b = indx + [b[-1]]   # two stage [vert1, vert2, [lamc1,lamc2]]
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

                dterms = [da1d12,db1d12,dg1d12,da2d12,db2d12,dg2d12]

                dtd12  = self.ambientSpace.geometry.funcDict["dtd12"](dataList.data[b[0]],dataList.data[b[1]],dterms)

                # a1,b1,g1 = dataList.data[b[0]].pos.copy()
                # a2,b2,g2 = dataList.data[b[1]].pos.copy()
                # print(type(b[-1]))
                riga1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[0,0], b[-1][stage], d12,  da1d12)
                rigb1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[1,1], b[-1][stage], d12,  db1d12)
                rigg1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[2,2], b[-1][stage], d12,  dg1d12)
                riga2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[0,0], b[-1][stage], d12,  da2d12)
                rigb2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[1,1], b[-1][stage], d12,  db2d12)
                rigg2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[2,2], b[-1][stage], d12,  dg2d12)


                rig_con   = self.ambientSpace.geometry.funcDict["rig_con"](self.x, d12)
                dtrig_con = self.ambientSpace.geometry.funcDict["rig_con_tderivative1"](self.x, d12, dtd12)

                rigdStates = [dState(zeros(3),array([riga1,rigb1,rigg1])), dState(zeros(3),array([riga2,rigb2,rigg2]))]

                temparr[b[0]] = temparr[b[0]].sub(rigdStates[0])
                temparr[b[1]] = temparr[b[1]].sub(rigdStates[1])

                # print("rig_con")
                # print(rig_con)
                tempconarr.append(rig_con)
                tempconarr.append(dtrig_con)

        # res_array = []
        # for c in temparr:
        #     res_array.append(c.vel.copy())
        #     res_array.append(c.acc.copy())

        # res_array = np.array(res_array).flatten

        # return [DataList(temparr, connectivity=dataList.connectivity), res_array]
        return [DataList(temparr, connectivity=dataList.connectivity, rig_connectivity=dataList.rig_connectivity),np.array(tempconarr)]

    # Jacobian will be constructed to be largely diagonal to help with efficiently for larger (sparse matrix) systems - block diagonal jacobian
    # This is in contrast with previous iteration of the solver design where top half was velocity residual expressions and bottom half were acceleration residual expressions
    def dynjac(self, dataList, stage, kList, params = None):
        # Dimension of ambient space (hard coded to 3D at the moment)
        ambient_dim = 3
        # Number of vertices in system
        vert_num = len(dataList.data)
        # Number of constraints
        con_num = len(dataList.rig_connectivity)*2 # Position and velocity constraints
        # Initialize jacobian matrix (multiply by 2 since considering residuals of velocity and accelerations for each vertex)
        jacobian = np.zeros((2*ambient_dim*vert_num + con_num,2*ambient_dim*vert_num + con_num))

        # Populate vertex data (block diagonal information in jacobian)
        for a in range(vert_num):
            tempmat = self.ambientSpace.geometry.funcDict["vertex_jac"](dataList.data[a])
            # print(a*2*ambient_dim )
            # print(a*2*ambient_dim + 2 * ambient_dim)
            jacobian[a*2*ambient_dim : a*2*ambient_dim + 2 * ambient_dim , a*2*ambient_dim : a*2*ambient_dim + 2 * ambient_dim] = tempmat.copy()

        # ks,x = [1,1]
        # m = 1.

        # Contributions from coupling potential (right now only pairwise coupling)
        if len(dataList.rig_connectivity) != 0:
            for b in dataList.rig_connectivity:
                # Change format to be in the form of 2-stage radau method
                if type(b[-1]).__name__ != 'list':
                    # print("reformat")
                    b[-1] = [b[-1],b[-1]]
                # print(b)
                indx = b[0:2]
                # print(indx)
                indx.sort()
                b = indx + [b[-1]]   # two stage [vert1, vert2, [lamc1,lamc2]]
                # b.sort()
                # Distance Function
                d12 = self.ambientSpace.geometry.funcDict["d12"](dataList.data[b[0]],dataList.data[b[1]])
                # First derivatives of distance function
                da1d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[0]],dataList.data[b[1]])
                db1d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[0]],dataList.data[b[1]])
                dg1d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[0]],dataList.data[b[1]])
                da2d12 = self.ambientSpace.geometry.funcDict["da1d12"](dataList.data[b[1]],dataList.data[b[0]])
                db2d12 = self.ambientSpace.geometry.funcDict["db1d12"](dataList.data[b[1]],dataList.data[b[0]])
                dg2d12 = self.ambientSpace.geometry.funcDict["dg1d12"](dataList.data[b[1]],dataList.data[b[0]])

                dterms = [da1d12,db1d12,dg1d12,da2d12,db2d12,dg2d12]

                if stage==1:
                    dtd12  = self.ambientSpace.geometry.funcDict["dtd12" ](kList.data[b[0]],kList.data[b[1]],dterms)
                else:
                    dtd12 = 0.

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

                ddterms = [
                    [da1d12a1,db1d12a1,dg1d12a1, da2d12a1,db2d12a1,dg2d12a1],
                    [da1d12b1,db1d12b1,dg1d12b1, da2d12b1,db2d12b1,dg2d12b1],
                    [da1d12g1,db1d12g1,dg1d12g1, da2d12g1,db2d12g1,dg2d12g1],

                    [da2d12a2,db1d12a2,dg1d12a2, da2d12a2,db2d12a2,dg2d12a2],
                    [da2d12b2,db1d12b2,dg1d12b2, da2d12b2,db2d12b2,dg2d12b2],
                    [da2d12g2,db1d12g2,dg1d12g2, da2d12g2,db2d12g2,dg2d12g2]
                    ]

                # if stage ==1:
                #      ddtd12  = self.ambientSpace.geometry.funcDict["ddtd12" ](kList.data[b[0]],kList.data[b[1]],dterms,ddterms)
                # else:
                #     ddtd12 = 0.

                # lagrange multplier derivatives
                dlamriga1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[0,0], 1., d12,  da1d12)
                dlamrigb1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[1,1], 1., d12,  db1d12)
                dlamrigg1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[2,2], 1., d12,  dg1d12)
                dlamriga2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[0,0], 1., d12,  da2d12)
                dlamrigb2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[1,1], 1., d12,  db2d12)
                dlamrigg2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[2,2], 1., d12,  dg2d12)

                metric1diag = self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos.copy()).diagonal()
                metric2diag = self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos.copy()).diagonal()

                dmetric_terms1 = self.ambientSpace.geometry.funcDict["dmetric_terms"](dataList.data[b[0]])
                dmetric_terms2 = self.ambientSpace.geometry.funcDict["dmetric_terms"](dataList.data[b[1]])

                # a1,b1,g1 = dataList.data[b[0]].pos.copy()
                # a2,b2,g2 = dataList.data[b[1]].pos.copy()

                #---------- particle 1

                da1da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1][stage], d12, da1d12, da1d12, dmetric_terms1[0,0], da1d12a1)
                db1da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1][stage], d12, da1d12, db1d12, dmetric_terms1[0,1], db1d12a1)
                dg1da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1][stage], d12, da1d12, dg1d12, dmetric_terms1[0,2], dg1d12a1)

                da2da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1][stage], d12, da1d12, da2d12, dmetric_terms1[0,3], da2d12a1)
                db2da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1][stage], d12, da1d12, db2d12, dmetric_terms1[0,4], db2d12a1)
                dg2da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1][stage], d12, da1d12, dg2d12, dmetric_terms1[0,5], dg2d12a1)

                #----------

                da1db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1][stage], d12, db1d12, da1d12, dmetric_terms1[1,0], da1d12b1)
                db1db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1][stage], d12, db1d12, db1d12, dmetric_terms1[1,1], db1d12b1)
                dg1db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1][stage], d12, db1d12, dg1d12, dmetric_terms1[1,2], dg1d12b1)

                da2db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1][stage], d12, db1d12, da2d12, dmetric_terms1[1,3], da2d12b1)
                db2db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1][stage], d12, db1d12, db2d12, dmetric_terms1[1,4], db2d12b1)
                dg2db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1][stage], d12, db1d12, dg2d12, dmetric_terms1[1,5], dg2d12b1)

                #------------

                da1dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1][stage], d12, dg1d12, da1d12, dmetric_terms1[2,0], da1d12g1)
                db1dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1][stage], d12, dg1d12, db1d12, dmetric_terms1[2,1], db1d12g1)
                dg1dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1][stage], d12, dg1d12, dg1d12, dmetric_terms1[2,2], dg1d12g1)

                da2dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1][stage], d12, dg1d12, da2d12, dmetric_terms1[2,3], da2d12g1)
                db2dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1][stage], d12, dg1d12, db2d12, dmetric_terms1[2,4], db2d12g1)
                dg2dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1][stage], d12, dg1d12, dg2d12, dmetric_terms1[2,5], dg2d12g1)


                # Contribution to particle 1 dynamics (block diagonal contribution)
                part1_geoterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da1da1con12,db1da1con12,dg1da1con12,  0., 0., 0.],
                    [da1db1con12,db1db1con12,dg1db1con12,  0., 0., 0.],
                    [da1dg1con12,db1dg1con12,dg1dg1con12,  0., 0., 0.]
                ])

                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] + part1_geoterms.copy()

                # # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                # part1_spterms = np.array([
                #     [da2da1con12,db2da1con12,dg2da1con12,  0., 0., 0.],
                #     [da2db1con12,db2db1con12,dg2db1con12,  0., 0., 0.],
                #     [da2dg1con12,db2dg1con12,dg2dg1con12,  0., 0., 0.],
                #     [0., 0., 0.,                           0., 0., 0.],
                #     [0., 0., 0.,                           0., 0., 0.],
                #     [0., 0., 0.,                           0., 0., 0.]
                # ])

                # jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part1_spterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part1_conterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da2da1con12,db2da1con12,dg2da1con12,  0., 0., 0.],
                    [da2db1con12,db2db1con12,dg2db1con12,  0., 0., 0.],
                    [da2dg1con12,db2dg1con12,dg2dg1con12,  0., 0., 0.]
                ])

                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part1_conterms.copy()

                # Constraint contributions (lambda columns)
                part1_langcolterms = np.array([0.,0.,0.,dlamriga1,dlamrigb1,dlamrigg1])

                # print(jacobian.shape)
                # print(jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + (stage)].shape)
                # print(part1_langcolterms)    
                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + (stage)] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + (stage)] + part1_langcolterms.copy()

                #------------- particle 2

                da1da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1][stage], d12, da2d12, da1d12, dmetric_terms2[0,0], da1d12a2)
                db1da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1][stage], d12, da2d12, db1d12, dmetric_terms2[0,1], db1d12a2)
                dg1da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1][stage], d12, da2d12, dg1d12, dmetric_terms2[0,2], dg1d12a2)

                da2da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1][stage], d12, da2d12, da2d12, dmetric_terms2[0,3], da2d12a2)
                db2da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1][stage], d12, da2d12, db2d12, dmetric_terms2[0,4], db2d12a2)
                dg2da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1][stage], d12, da2d12, dg2d12, dmetric_terms2[0,5], dg2d12a2)

                #----------

                da1db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1][stage], d12, db2d12, da1d12, dmetric_terms2[1,0], da1d12b2)
                db1db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1][stage], d12, db2d12, db1d12, dmetric_terms2[1,1], db1d12b2)
                dg1db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1][stage], d12, db2d12, dg1d12, dmetric_terms2[1,2], dg1d12b2)

                da2db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1][stage], d12, db2d12, da2d12, dmetric_terms2[1,3], da2d12b2)
                db2db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1][stage], d12, db2d12, db2d12, dmetric_terms2[1,4], db2d12b2)
                dg2db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1][stage], d12, db2d12, dg2d12, dmetric_terms2[1,5], dg2d12b2)

                #------------

                da1dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1][stage], d12, dg2d12, da1d12, dmetric_terms2[2,0], da1d12g2)
                db1dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1][stage], d12, dg2d12, db1d12, dmetric_terms2[2,1], db1d12g2)
                dg1dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1][stage], d12, dg2d12, dg1d12, dmetric_terms2[2,2], dg1d12g2)

                da2dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1][stage], d12, dg2d12, da2d12, dmetric_terms2[2,3], da2d12g2)
                db2dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1][stage], d12, dg2d12, db2d12, dmetric_terms2[2,4], db2d12g2)
                dg2dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1][stage], d12, dg2d12, dg2d12, dmetric_terms2[2,5], dg2d12g2)


                # Contribution to particle 1 dynamics (block diagonal contribution)
                part2_geoterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da2da2con12,db2da2con12,dg2da2con12,  0., 0., 0.],
                    [da2db2con12,db2db2con12,dg2db2con12,  0., 0., 0.],
                    [da2dg2con12,db2dg2con12,dg2dg2con12,  0., 0., 0.]
                ])

                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part2_geoterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part2_spterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da1da2con12,db1da2con12,dg1da2con12,  0., 0., 0.],
                    [da1db2con12,db1db2con12,dg1db2con12,  0., 0., 0.],
                    [da1dg2con12,db1dg2con12,dg1dg2con12,  0., 0., 0.]
                ])

                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] + part2_spterms.copy()

                # Constraint contributions (lambda columns)
                part2_langcolterms = np.array([0.,0.,0.,dlamriga2,dlamrigb2,dlamrigg2])

                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + (stage)] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + (stage)] + part2_langcolterms.copy()

                # Constraint contributions (lambda rows) - only for stage s for s-stage method
                if stage == 1:
                    langrowconterms = np.array([dlamriga1,dlamrigb1,dlamrigg1, 0.,0.,0., dlamriga2,dlamrigb2,dlamrigg2, 0.,0.,0.])

                    da1langrowdtconterms = self.ambientSpace.geometry.funcDict["rig_dtcon_derivative1_array"](dataList.data[b[0]],dataList.data[b[1]],ddterms)
                    dda1langrowdtconterms = np.array([1.,1.,1., 1.,1.,1.]) @ ddterms

                    langrowdtconterms = np.zeros(2*vert_num*ambient_dim)

                    langrowdtconterms[0:3] = da1langrowdtconterms[0:3]
                    langrowdtconterms[3:6] = dda1langrowdtconterms[0:3]
                    langrowdtconterms[0:3] = da1langrowdtconterms[3:6]
                    langrowdtconterms[3:6] = dda1langrowdtconterms[3:6]

                    jacobian[vert_num*2*ambient_dim , 0 : vert_num*2*ambient_dim] = jacobian[vert_num*2*ambient_dim , 0 : vert_num*2*ambient_dim] +  langrowconterms.copy()
                    jacobian[vert_num*2*ambient_dim + 1 , 0 : vert_num*2*ambient_dim] = jacobian[vert_num*2*ambient_dim + 1 , 0 : vert_num*2*ambient_dim] +  langrowdtconterms.copy()


        # fullj = np.eye(2*ambient_dim*vert_num) - jacobian
        return [jacobian[0 : vert_num*2*ambient_dim , 0 : vert_num*2*ambient_dim], jacobian]

    def difffunc(self, dataList, k1, k2, conterms):
        # Set to run Gauss 2-stage method
        a11,a12 = [5./12., -1./12.]
        a21,a22 = [3./4., 1./4.]

        return np.array([
            (k1.toArray() - self.dynfunc(self.arrayToDataList(dataList.toArray() + (a11*k1.toArray() + a12*k2.toArray())*self.stepSize, dataList), stage=0)[0].toArray()).tolist()+
            (k2.toArray() - self.dynfunc(self.arrayToDataList(dataList.toArray() + (a21*k1.toArray() + a22*k2.toArray())*self.stepSize, dataList), stage=1)[0].toArray()).tolist()+
            (conterms).tolist()
        ]).flatten()

    def step(self, dataList, tol = 1e-15, imax = 10):
        # Dimension of ambient space (hard coded to 3D at the moment)
        ambient_dim = 3
        # Number of vertices in system
        vert_num = len(dataList.data)
        # Number of constraints
        con_num = len(dataList.rig_connectivity)*2 # Position and velocity constraints

        # Set to run Gauss 2-stage method
        a11,a12 = [5./12., -1./12.]
        a21,a22 = [3./4., 1./4.]
        bs1,bs2 = [3./4., 1./4.]

        # Initial Guess - Explicit Euler
        # k, conterms = self.dynfunc(dataList, stage=1)
        # x1guess = dataList.toArray() + (1./3.)*self.stepSize*k.toArray()
        # x2guess = dataList.toArray() + (1.)*self.stepSize*k.toArray()
        # k1, _        = self.dynfunc(self.arrayToDataList(x1guess, dataList), stage=0)
        # k2, conterms = self.dynfunc(self.arrayToDataList(x2guess, dataList), stage=1) # constraints only satisfied at last stage
        k1,   _       = self.dynfunc(dataList, stage=0)
        k2, conterms2 = self.dynfunc(dataList, stage=1) # constraints only satisfied at last stage
        print("velocity")
        print(k1.data[0].vel)
        print(k1.data[1].vel)
        print("acceleration")
        print(k1.data[0].acc)
        print(k1.data[1].acc)
        print(conterms2)

        # Check Error Before iterations
        er = self.difffunc(dataList, k1, k2, conterms2)
        print(er)

        # Begin Iterations
        for a in range(imax):
            if np.linalg.norm(er) >= tol:
                j1, j1wcon = self.dynjac(self.arrayToDataList(dataList.toArray() + (a11*k1.toArray() + a12*k2.toArray())*self.stepSize, dataList), stage=0, kList = k1)
                j2, j2wcon = self.dynjac(self.arrayToDataList(dataList.toArray() + (a21*k1.toArray() + a22*k2.toArray())*self.stepSize, dataList), stage=1, kList = k2)
                
                fullj = np.block([
                    [np.eye(j1.shape[0]) - self.stepSize*a11*j1, -self.stepSize*a12*j1,-self.stepSize*a12*j1wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num]],
                    [-self.stepSize*a21*j2, np.eye(j2.shape[0]) - self.stepSize*a22*j2,-self.stepSize*a12*j2wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num]],
                    [np.zeros((con_num,j1.shape[1])),j2wcon[vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num , 0 : vert_num*2*ambient_dim],np.zeros((con_num,con_num))]
                ])

                print(fullj)
                linsolve = np.linalg.solve(fullj,-er)
                print(linsolve)

                k1 = self.arrayToDataList(k1.toArray() + linsolve[0:j1.shape[0]],dataList)
                k2 = self.arrayToDataList(k2.toArray() + linsolve[j2.shape[0]:2*j2.shape[0]],dataList)
                # print(linsolve[2*j2.shape[0]:2*j2.shape[0]+con_num])
                for b in dataList.rig_connectivity:
                    b[-1] = (np.array(b[-1]) + linsolve[2*j2.shape[0]:2*j2.shape[0]+con_num]).tolist()
                print(dataList.rig_connectivity)

                er = self.difffunc(dataList, k1, k2, conterms2)
            else:
                break


        newarr = dataList.toArray() + self.stepSize*(bs1*k1.toArray() + bs2*k2.toArray())
        return self.arrayToDataList(newarr, dataList)