# // a class for integrating equations of motion called "integrator" and one specific implementation, Runge Kutta
# //derive is a function taking a state to state (now storing velocity and acceleration instead of position and velocity)
# //items fed into RungeKutta need to have the following methods available:
# //.add, .multiplyScalar, .clone

#######################
## Currently not functional!!!!!!
#######################

import numpy as np
from copy import deepcopy
from numpy import sin, cos, tan, cosh, sinh, tanh, arccosh, sqrt, array, zeros

from src.Computation.DataList import DataList
from src.Computation.State import State
from src.Computation.dState import dState

# //implementing the Rk4 Scheme for arbitrary classes that have clone add and multiplyScalar
# //will use this possibly on individual states, or on entire DataLists!
class RigidRadau3:
    """
    A class used to perform numerical integration via Implicit 3-stage Radau method with rigidity constraints(RigidRS3)

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
        chunkarr = []
        # print(arr)
        # Separate into chucks of 3 values
        for i in range(0, len(arr), 3):
            chunkarr.append(arr[i:i + 3])
        # print(chunkarr)
            
        # For DataList of States
        if dataList.data[0].__class__.__name__ == "State":
            for a in range(0,len(chunkarr),2):
                # print(chunkarr[a])
                # print(chunkarr[a+1])
                temparr.append(State(chunkarr[a],chunkarr[a+1]))

        # For DataList of dStates
        if dataList.data[0].__class__.__name__ == "dState":
            for a in range(0,len(chunkarr),2):
                # print(chunkarr[a])
                # print(chunkarr[a+1])
                temparr.append(dState(chunkarr[a],chunkarr[a+1]))

        # for a in range(0,len(chunkarr),2):
        #     # print(chunkarr[a])
        #     # print(chunkarr[a+1])
        #     temparr.append(State(chunkarr[a],chunkarr[a+1]))

        # for a in range(len(dataList.data)):
        #     # print(a)
        #     temparr.append(State(arr[a*2*3:a*2*3 + 3],arr[a*2*3 + 3:a*2*3 + 6]))

        return DataList(temparr, connectivity=deepcopy(dataList.connectivity), rig_connectivity=deepcopy(dataList.rig_connectivity))

    def dynfunc(self, dataList, params = None):
        temparr = []
        tempconarr = []
        tempdtconarr = []
        tempddtconarr = []
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
            # First pass for position and velocity constraint
            for b in dataList.rig_connectivity:
                indx = b[0:2]
                indx.sort()
                b = indx + [b[-1]]   # two stage [vert1, vert2, lam]
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

                
                dtd12  = self.ambientSpace.geometry.funcDict["dtd12" ](dataList.data[b[0]],dataList.data[b[1]],dterms)


                # a1,b1,g1 = dataList.data[b[0]].pos.copy()
                # a2,b2,g2 = dataList.data[b[1]].pos.copy()
                # print(type(b[-1]))
                riga1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[0,0], b[-1], d12,  da1d12)
                rigb1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[1,1], b[-1], d12,  db1d12)
                rigg1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[2,2], b[-1], d12,  dg1d12)
                riga2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[0,0], b[-1], d12,  da2d12)
                rigb2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[1,1], b[-1], d12,  db2d12)
                rigg2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[2,2], b[-1], d12,  dg2d12)

                # print([riga1,rigb1,rigg1])

                rigdStates = [dState(zeros(3),array([riga1,rigb1,rigg1])), dState(zeros(3),array([riga2,rigb2,rigg2]))]

                temparr[b[0]] = temparr[b[0]].sub(rigdStates[0])
                temparr[b[1]] = temparr[b[1]].sub(rigdStates[1])
                # print(temparr[b[0]])
                # print(temparr[b[1]])

                rig_con    = self.ambientSpace.geometry.funcDict["rig_con"](self.x, d12)
                dtrig_con  = self.ambientSpace.geometry.funcDict["rig_con_tderivative1"](self.x, d12, dtd12)
                # ddtrig_con = self.ambientSpace.geometry.funcDict["rig_con_tderivative2"](self.x, d12, dtd12, ddtd12)

                # print("rig_con")
                # print(rig_con)
                tempconarr.append(rig_con)
                tempdtconarr.append(dtrig_con)
                # tempddtconarr.append(ddtrig_con)

            # Second pass for acceleration constraint (Needs optimization)
            for b in dataList.rig_connectivity:
                indx = b[0:2]
                indx.sort()
                b = indx + [b[-1]]   # two stage [vert1, vert2, lam]
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

                
                dtd12  = self.ambientSpace.geometry.funcDict["dtd12" ](dataList.data[b[0]],dataList.data[b[1]],dterms)

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

                    [da1d12a2,db1d12a2,dg1d12a2, da2d12a2,db2d12a2,dg2d12a2],
                    [da1d12b2,db1d12b2,dg1d12b2, da2d12b2,db2d12b2,dg2d12b2],
                    [da1d12g2,db1d12g2,dg1d12g2, da2d12g2,db2d12g2,dg2d12g2]
                    ]


                # temparr[b[0]] = temparr[b[0]].sub(rigdStates[0])
                # temparr[b[1]] = temparr[b[1]].sub(rigdStates[1])
                # print(temparr[b[0]])
                # print(temparr[b[1]])

                ddtd12  = self.ambientSpace.geometry.funcDict["ddtd12" ](temparr[b[0]],temparr[b[1]],dterms,ddterms)

                # rig_con    = self.ambientSpace.geometry.funcDict["rig_con"](self.x, d12)
                # dtrig_con  = self.ambientSpace.geometry.funcDict["rig_con_tderivative1"](self.x, d12, dtd12)
                ddtrig_con = self.ambientSpace.geometry.funcDict["rig_con_tderivative2"](self.x, d12, dtd12, ddtd12)

                # print("rig_con")
                # print(rig_con)
                # tempconarr.append(rig_con)
                # tempdtconarr.append(dtrig_con)
                tempddtconarr.append(ddtrig_con)
            

        # res_array = []
        # for c in temparr:
        #     res_array.append(c.vel.copy())
        #     res_array.append(c.acc.copy())

        # res_array = np.array(res_array).flatten

        # return [DataList(temparr, connectivity=dataList.connectivity), res_array]
        return [DataList(temparr, connectivity=deepcopy(dataList.connectivity), rig_connectivity=deepcopy(dataList.rig_connectivity)),np.array(tempconarr),np.array(tempdtconarr),np.array(tempddtconarr)]

    # Jacobian will be constructed to be largely diagonal to help with efficiently for larger (sparse matrix) systems - block diagonal jacobian
    # This is in contrast with previous iteration of the solver design where top half was velocity residual expressions and bottom half were acceleration residual expressions
    def dynjac(self, dataList, params = None):
        # Dimension of ambient space (hard coded to 3D at the moment)
        ambient_dim = 3
        # Number of vertices in system
        vert_num = len(dataList.data)
        # Number of constraints
        con_num = int(len(dataList.rig_connectivity)*2) # Position and velocity constraints both share lagrange multiplier
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
                indx = b[0:2]
                # print(indx)
                indx.sort()
                b = indx + [b[-1]]   # two stage [vert1, vert2, lam]
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

                
                dtd12  = self.ambientSpace.geometry.funcDict["dtd12" ](dataList.data[b[0]],dataList.data[b[1]],dterms)

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

                # lagrange multplier derivatives for columns
                dlamriga1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[0,0], 1., d12,  da1d12)
                dlamrigb1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[1,1], 1., d12,  db1d12)
                dlamrigg1 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos)[2,2], 1., d12,  dg1d12)
                dlamriga2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[0,0], 1., d12,  da2d12)
                dlamrigb2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[1,1], 1., d12,  db2d12)
                dlamrigg2 = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](self.m, self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos)[2,2], 1., d12,  dg2d12)

                # lagrange multplier derivatives for position constraint row
                da1lamcon = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](1., 1., 1., d12,  da1d12)
                db1lamcon = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](1., 1., 1., d12,  db1d12)
                dg1lamcon = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](1., 1., 1., d12,  dg1d12)
                da2lamcon = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](1., 1., 1., d12,  da2d12)
                db2lamcon = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](1., 1., 1., d12,  db2d12)
                dg2lamcon = self.ambientSpace.geometry.funcDict["rig_con_derivative1"](1., 1., 1., d12,  dg2d12)

                # lagrange multplier derivatives for velocity constraint row
                da1langdtconarr = self.ambientSpace.geometry.funcDict["rig_dtcon_derivative1_array"](dataList.data[b[0]],dataList.data[b[1]], d12, dtd12, dterms, ddterms)

                metric1diag = self.ambientSpace.geometry.metricTensor(dataList.data[b[0]].pos.copy()).diagonal()
                metric2diag = self.ambientSpace.geometry.metricTensor(dataList.data[b[1]].pos.copy()).diagonal()

                dmetric_terms1 = self.ambientSpace.geometry.funcDict["dmetric_terms"](dataList.data[b[0]])
                dmetric_terms2 = self.ambientSpace.geometry.funcDict["dmetric_terms"](dataList.data[b[1]])

                # a1,b1,g1 = dataList.data[b[0]].pos.copy()
                # a2,b2,g2 = dataList.data[b[1]].pos.copy()

                #---------- particle 1

                da1da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1], d12, da1d12, da1d12, dmetric_terms1[0,0], da1d12a1)
                db1da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1], d12, da1d12, db1d12, dmetric_terms1[0,1], db1d12a1)
                dg1da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1], d12, da1d12, dg1d12, dmetric_terms1[0,2], dg1d12a1)

                da2da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1], d12, da1d12, da2d12, dmetric_terms1[0,3], da2d12a1)
                db2da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1], d12, da1d12, db2d12, dmetric_terms1[0,4], db2d12a1)
                dg2da1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[0], b[-1], d12, da1d12, dg2d12, dmetric_terms1[0,5], dg2d12a1)

                #----------

                da1db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1], d12, db1d12, da1d12, dmetric_terms1[1,0], da1d12b1)
                db1db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1], d12, db1d12, db1d12, dmetric_terms1[1,1], db1d12b1)
                dg1db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1], d12, db1d12, dg1d12, dmetric_terms1[1,2], dg1d12b1)

                da2db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1], d12, db1d12, da2d12, dmetric_terms1[1,3], da2d12b1)
                db2db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1], d12, db1d12, db2d12, dmetric_terms1[1,4], db2d12b1)
                dg2db1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[1], b[-1], d12, db1d12, dg2d12, dmetric_terms1[1,5], dg2d12b1)

                #------------

                da1dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1], d12, dg1d12, da1d12, dmetric_terms1[2,0], da1d12g1)
                db1dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1], d12, dg1d12, db1d12, dmetric_terms1[2,1], db1d12g1)
                dg1dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1], d12, dg1d12, dg1d12, dmetric_terms1[2,2], dg1d12g1)

                da2dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1], d12, dg1d12, da2d12, dmetric_terms1[2,3], da2d12g1)
                db2dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1], d12, dg1d12, db2d12, dmetric_terms1[2,4], db2d12g1)
                dg2dg1con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric1diag[2], b[-1], d12, dg1d12, dg2d12, dmetric_terms1[2,5], dg2d12g1)


                #Contribution to particle 1 dynamics (block diagonal contribution)
                part1_geoterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da1da1con12,db1da1con12,dg1da1con12,  0., 0., 0.],
                    [da1db1con12,db1db1con12,dg1db1con12,  0., 0., 0.],
                    [da1dg1con12,db1dg1con12,dg1dg1con12,  0., 0., 0.]
                ])

                # Minus for correct sign when considering residuals
                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] - part1_geoterms.copy()

                # # # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                # # part1_spterms = np.array([
                # #     [da2da1con12,db2da1con12,dg2da1con12,  0., 0., 0.],
                # #     [da2db1con12,db2db1con12,dg2db1con12,  0., 0., 0.],
                # #     [da2dg1con12,db2dg1con12,dg2dg1con12,  0., 0., 0.],
                # #     [0., 0., 0.,                           0., 0., 0.],
                # #     [0., 0., 0.,                           0., 0., 0.],
                # #     [0., 0., 0.,                           0., 0., 0.]
                # # ])

                # # jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] + part1_spterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part1_conterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da2da1con12,db2da1con12,dg2da1con12,  0., 0., 0.],
                    [da2db1con12,db2db1con12,dg2db1con12,  0., 0., 0.],
                    [da2dg1con12,db2dg1con12,dg2dg1con12,  0., 0., 0.]
                ])

                # Minus for correct sign when considering residuals
                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] - part1_conterms.copy()

                # Constraint contributions (lambda columns)
                part1_langcolterms = np.array([0.,0.,0.,dlamriga1,dlamrigb1,dlamrigg1])
                # print("langcolterms")
                # print(part1_langcolterms)

                # print(jacobian.shape)
                # print(jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + (stage)].shape)
                # print(part1_langcolterms)    
                # Minus for correct sign when considering residuals
                jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + dataList.rig_connectivity.index(b)] = jacobian[b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + dataList.rig_connectivity.index(b)] - part1_langcolterms.copy()

                #------------- particle 2

                da1da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1], d12, da2d12, da1d12, dmetric_terms2[0,0], da1d12a2)
                db1da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1], d12, da2d12, db1d12, dmetric_terms2[0,1], db1d12a2)
                dg1da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1], d12, da2d12, dg1d12, dmetric_terms2[0,2], dg1d12a2)

                da2da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1], d12, da2d12, da2d12, dmetric_terms2[0,3], da2d12a2)
                db2da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1], d12, da2d12, db2d12, dmetric_terms2[0,4], db2d12a2)
                dg2da2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[0], b[-1], d12, da2d12, dg2d12, dmetric_terms2[0,5], dg2d12a2)

                #----------

                da1db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1], d12, db2d12, da1d12, dmetric_terms2[1,0], da1d12b2)
                db1db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1], d12, db2d12, db1d12, dmetric_terms2[1,1], db1d12b2)
                dg1db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1], d12, db2d12, dg1d12, dmetric_terms2[1,2], dg1d12b2)

                da2db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1], d12, db2d12, da2d12, dmetric_terms2[1,3], da2d12b2)
                db2db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1], d12, db2d12, db2d12, dmetric_terms2[1,4], db2d12b2)
                dg2db2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[1], b[-1], d12, db2d12, dg2d12, dmetric_terms2[1,5], dg2d12b2)

                #------------

                da1dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1], d12, dg2d12, da1d12, dmetric_terms2[2,0], da1d12g2)
                db1dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1], d12, dg2d12, db1d12, dmetric_terms2[2,1], db1d12g2)
                dg1dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1], d12, dg2d12, dg1d12, dmetric_terms2[2,2], dg1d12g2)

                da2dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1], d12, dg2d12, da2d12, dmetric_terms2[2,3], da2d12g2)
                db2dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1], d12, dg2d12, db2d12, dmetric_terms2[2,4], db2d12g2)
                dg2dg2con12 = self.ambientSpace.geometry.funcDict["rig_con_derivative2"](self.m, metric2diag[2], b[-1], d12, dg2d12, dg2d12, dmetric_terms2[2,5], dg2d12g2)


                #Contribution to particle 1 dynamics (block diagonal contribution)
                part2_geoterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da2da2con12,db2da2con12,dg2da2con12,  0., 0., 0.],
                    [da2db2con12,db2db2con12,dg2db2con12,  0., 0., 0.],
                    [da2dg2con12,db2dg2con12,dg2dg2con12,  0., 0., 0.]
                ])

                # Minus for correct sign when considering residuals
                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] - part2_geoterms.copy()

                # # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                # part2_spterms = np.array([
                #     [0., 0., 0.,                           0., 0., 0.],
                #     [0., 0., 0.,                           0., 0., 0.],
                #     [0., 0., 0.,                           0., 0., 0.],
                #     [da1da2con12,db1da2con12,dg1da2con12,  0., 0., 0.],
                #     [da1db2con12,db1db2con12,dg1db2con12,  0., 0., 0.],
                #     [da1dg2con12,db1dg2con12,dg1dg2con12,  0., 0., 0.]
                # ])

                # # Minus for correct sign when considering residuals
                # jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] - part2_spterms.copy()

                # Contribution to particle 1 dynamics from particle 2 (off diagonal contribution)
                part2_conterms = np.array([
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [0., 0., 0.,                           0., 0., 0.],
                    [da1da2con12,db1da2con12,dg1da2con12,  0., 0., 0.],
                    [da1db2con12,db1db2con12,dg1db2con12,  0., 0., 0.],
                    [da1dg2con12,db1dg2con12,dg1dg2con12,  0., 0., 0.]
                ])

                # Minus for correct sign when considering residuals
                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] - part2_conterms.copy()

                # Constraint contributions (lambda columns)
                part2_langcolterms = np.array([0.,0.,0.,dlamriga2,dlamrigb2,dlamrigg2])
                # print(part2_langcolterms)

                # Minus for correct sign when considering residuals
                jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + dataList.rig_connectivity.index(b)] = jacobian[b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim , vert_num*2*ambient_dim + dataList.rig_connectivity.index(b)] - part2_langcolterms.copy()

                # Constraint contributions (lambda rows) - only for stage s for s-stage method
                langrowconterms = np.array([da1lamcon,db1lamcon,dg1lamcon, 0.,0.,0., da2lamcon,db2lamcon,dg2lamcon, 0.,0.,0.])
                # print("gamma term")
                # print(dg1lamcon)

                # da1langrowdtconterms = self.ambientSpace.geometry.funcDict["rig_dtcon_derivative1_array"](dataList.data[b[0]],dataList.data[b[1]],ddterms)
                # dda1langrowdtconterms = np.array([1.,1.,1., 1.,1.,1.]) @ ddterms

                langrowdtconterms = da1langdtconarr
                # print(b)
                # print(langrowconterms)
                # print(langrowdtconterms)

                # langrowdtconterms[0:3] = da1langrowdtconterms[0:3]
                # langrowdtconterms[3:6] = dda1langrowdtconterms[0:3]
                # langrowdtconterms[0:3] = da1langrowdtconterms[3:6]
                # langrowdtconterms[3:6] = dda1langrowdtconterms[3:6]

                jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b), b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b), b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] +  langrowconterms.copy()[0:6]
                jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b), b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b), b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] +  langrowconterms.copy()[6:12]

                jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b) + con_num // 2, b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] = jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b) + con_num // 2, b[0]*2*ambient_dim : b[0]*2*ambient_dim + 2 * ambient_dim] +  langrowdtconterms.copy()[0:6]
                jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b) + con_num // 2, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] = jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b) + con_num // 2, b[1]*2*ambient_dim : b[1]*2*ambient_dim + 2 * ambient_dim] +  langrowdtconterms.copy()[6:12]

                # jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b), 0 : vert_num*2*ambient_dim] = jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b), 0 : vert_num*2*ambient_dim] +  langrowconterms.copy()
                # jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b) + con_num // 2 , 0 : vert_num*2*ambient_dim] = jacobian[vert_num*2*ambient_dim + dataList.rig_connectivity.index(b) + con_num // 2 , 0 : vert_num*2*ambient_dim] +  langrowdtconterms.copy()


        # fullj = np.eye(2*ambient_dim*vert_num) - jacobian
        return [jacobian[0 : vert_num*2*ambient_dim , 0 : vert_num*2*ambient_dim], jacobian]

    def difffunc(self, dataList, dataListk1n1, dataListk2n1, dataListk3n1, k1, k2, k3, conterms, dtconterms, ddtconterms):
        # Set to run Gauss 3-stage method
        root6 = np.sqrt(6.)
        # Set to run RadauIIA 3-stage method
        a11,a12,a13 = [(88.-7.*root6)/360., (296.-169.*root6)/1800., (-2.+3.*root6)/225.]
        a21,a22,a23 = [(296.+169.*root6)/1800., (88.+7.*root6)/360., (-2.-3.*root6)/225.]
        a31,a32,a33 = [(16.-root6)/36., (16.+root6)/36., 1./9.]

        # print(k1.toArray())
        

        return np.array([
            (dataListk1n1.toArray() - dataList.toArray() - self.stepSize*(a11*k1.toArray() + a12*k2.toArray() + a13*k3.toArray())).tolist()+
            (dataListk2n1.toArray() - dataList.toArray() - self.stepSize*(a21*k1.toArray() + a22*k2.toArray() + a23*k3.toArray())).tolist()+
            (dataListk3n1.toArray() - dataList.toArray() - self.stepSize*(a31*k1.toArray() + a32*k2.toArray() + a33*k3.toArray())).tolist()+
            (conterms).tolist()+
            (dtconterms).tolist()+
            (ddtconterms).tolist()
        ]).flatten()

    def step(self, dataList, tol = 1e-15, imax = 100):
        # Dimension of ambient space (hard coded to 3D at the moment)
        ambient_dim = 3
        # Number of vertices in system
        vert_num = int(len(dataList.data))
        # Number of constraints
        con_num = int(len(dataList.rig_connectivity)*3) # Position, velocity, and acceleration constraints

        # Set to run Gauss 3-stage method
        root6 = np.sqrt(6.)
        # Set to run RadauIIA 3-stage method
        a11,a12,a13 = [(88.-7.*root6)/360., (296.-169.*root6)/1800., (-2.+3.*root6)/225.]
        a21,a22,a23 = [(296.+169.*root6)/1800., (88.+7.*root6)/360., (-2.-3.*root6)/225.]
        a31,a32,a33 = [(16.-root6)/36., (16.+root6)/36., 1./9.]

        # Initial Guess - Explicit Euler
        k, _ , _ , _ = self.dynfunc(dataList)
        x1guess = dataList.toArray()
        x2guess = dataList.toArray()
        x3guess = dataList.toArray()
        x1guessdataList = self.arrayToDataList(x1guess, dataList)
        x2guessdataList = self.arrayToDataList(x2guess, dataList)
        x3guessdataList = self.arrayToDataList(x3guess, dataList)
        k1, _ , _  , _  = self.dynfunc(x1guessdataList.clone())
        k2, _ , _  , _  = self.dynfunc(x2guessdataList.clone())
        k3, conterms, dtconterms, ddtconterms = self.dynfunc(x3guessdataList.clone()) # constraints only satisfied at last stage



        # Check Error Before iterations
        er = self.difffunc(dataList, x1guessdataList, x2guessdataList, x3guessdataList, k1, k2, k3, conterms, dtconterms, ddtconterms)
        print("conlist")
        print(er)
        print("")

        # Begin Iterations
        for a in range(imax):
            if np.linalg.norm(er) >= tol:
                j1, j1wcon = self.dynjac(x1guessdataList.clone())
                j2, j2wcon = self.dynjac(x2guessdataList.clone())
                j3, j3wcon = self.dynjac(x3guessdataList.clone())

                # j1, j1wcon = self.dynjac(self.arrayToDataList(dataList.toArray() + (a11*k1.toArray() + a12*k2.toArray())*self.stepSize, dataList))
                # j2, j2wcon = self.dynjac(self.arrayToDataList(dataList.toArray() + (a21*k1.toArray() + a22*k2.toArray())*self.stepSize, dataList))

                lambda_block1 = np.zeros((j1.shape[0],con_num))
                lambda_block2 = np.zeros((j2.shape[0],con_num))
                # print(j1wcon.shape)
                # print(j1wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num//2] )
                # print(j1wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim + con_num//2 : vert_num*2*ambient_dim + con_num])

                lambda_block1[:,0      :con_num//2] = - self.stepSize*a11*j1wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num//2] 
                lambda_block1[:,con_num//2:con_num] = - self.stepSize*a12*j1wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num//2]
                lambda_block2[:,0      :con_num//2] = - self.stepSize*a21*j2wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num//2] 
                lambda_block2[:,con_num//2:con_num] = - self.stepSize*a22*j2wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num//2]
                
                # print(j1wcon[0 : vert_num*2*ambient_dim , vert_num*2*ambient_dim + con_num//2 : vert_num*2*ambient_dim + con_num].shape)
                # print(lambda_block2.shape)
                # print(lambda_block1)
                conblock = j2wcon[vert_num*2*ambient_dim : vert_num*2*ambient_dim + con_num , 0 : vert_num*2*ambient_dim]
                # print(conblock.shape)
                # print(conblock)
                
                fullj = np.block([
                    [np.eye(j1.shape[0]) - self.stepSize*a11*j1, - self.stepSize*a12*j1,lambda_block1],
                    [- self.stepSize*a21*j2, np.eye(j2.shape[0]) - self.stepSize*a22*j2,lambda_block2],
                    [np.zeros((con_num,j1.shape[1])),conblock,np.zeros((con_num,con_num))]
                ])

                # print("Full Jacobian")
                # print(fullj)
                # # fullj.astype('float64').tofile('np.dat')
                # print("")
                linsolve = np.linalg.solve(fullj,-er)
                # print("solution vector")
                # print(linsolve)

                x1guessdataList = self.arrayToDataList(x1guessdataList.toArray() + linsolve[0:j1.shape[0]],x1guessdataList)
                x2guessdataList = self.arrayToDataList(x2guessdataList.toArray() + linsolve[j2.shape[0]:2*j2.shape[0]],x2guessdataList)

                # print("rig_con")
                # print(x1guessdataList.rig_connectivity)
                # print(x2guessdataList.rig_connectivity)

                # print(linsolve[2*j2.shape[0] + con_num//2 : 2*j2.shape[0]+con_num])
                # print(linsolve[2*j2.shape[0]:2*j2.shape[0]+con_num])
                for b in range(len(dataList.rig_connectivity)):
                    # print(2*j2.shape[0] + b)
                    # print(linsolve[2*j2.shape[0] + b])
                    # print(np.array(x1guessdataList.rig_connectivity[b][-1]) + linsolve[2*j2.shape[0] + b])
                    x1guessdataList.rig_connectivity[b][-1] = (np.array(x1guessdataList.rig_connectivity[b][-1]) + linsolve[2*j2.shape[0] + (b)])
                    x2guessdataList.rig_connectivity[b][-1] = (np.array(x2guessdataList.rig_connectivity[b][-1]) + linsolve[2*j2.shape[0] + (b+con_num//2)])
                    # print(b)
                # print(dataList.rig_connectivity)
                # print(x1guessdataList.rig_connectivity)
                # print(x2guessdataList.rig_connectivity)

                k1, _ , _       = self.dynfunc(x1guessdataList.clone())
                k2, conterms, dtconterms = self.dynfunc(x2guessdataList.clone())  # constraints only satisfied at last stage

                # print("rig_con")
                # print(x1guessdataList.rig_connectivity)
                # print(x2guessdataList.rig_connectivity)



                # Need to update conterms

                er = self.difffunc(dataList, x1guessdataList, x2guessdataList, k1, k2, conterms, dtconterms)
                # print(a)
                # print("conlist")
                # print(er)
                # print("")
            else:
                break

        # print("rig_con")
        # print(x1guessdataList.rig_connectivity)
        # print(x2guessdataList.rig_connectivity)

        newarr = x2guessdataList.clone().toArray()
        return self.arrayToDataList(newarr, x2guessdataList)