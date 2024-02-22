# // a class for integrating equations of motion called "integrator" and one specific implementation, Runge Kutta
# //derive is a function taking a state to state (now storing velocity and acceleration instead of position and velocity)
# //items fed into RungeKutta need to have the following methods available:
# //.add, .multiplyScalar, .clone

from numpy import sin, cos, cosh, sinh, arccosh, sqrt, array, zeros

from src.Computation.DataList import DataList
from src.Computation.dState import dState

# //implementing the Rk4 Scheme for arbitrary classes that have clone add and multiplyScalar
# //will use this possibly on individual states, or on entire DataLists!
class RungeKutta:
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

    def deriveFunc(self, dataList):
        temparr = []

        # Contributions from ambient space geometry
        for a in range(len(dataList.data)):
            # print(a)
            temparr.append(self.ambientSpace.acceleration(dataList.data[a].clone()))

        ks,x = [1,1]
        m = 1.

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

                spa1 = (ks*(arccosh(d12) - x)*da1d12)/(m*sqrt(d12**2. - 1.))
                spb1 = (ks*(arccosh(d12) - x)*db1d12)/(m*sinh(a1)**2. * sqrt(d12**2. - 1.))
                spg1 = (ks*(arccosh(d12) - x)*dg1d12)/(m*sinh(a1)**2. * sin(b1)**2. * sqrt(d12**2. - 1.))
                spa2 = (ks*(arccosh(d12) - x)*da2d12)/(m*sqrt(d12**2. - 1.))
                spb2 = (ks*(arccosh(d12) - x)*db2d12)/(m*sinh(a2)**2. * sqrt(d12**2. - 1.))
                spg2 = (ks*(arccosh(d12) - x)*dg2d12)/(m*sinh(a2)**2. * sin(b2)**2. * sqrt(d12**2. - 1.))

                spdStates = [dState(zeros(3),array([spa1,spb1,spg1])), dState(zeros(3),array([spa2,spb2,spg2]))]

                # # sp1
                # print(spdStates[0].acc)
                # # sp2
                # print(spdStates[1].acc)
                temparr[b[0]] = temparr[b[0]].sub(spdStates[0])
                temparr[b[1]] = temparr[b[1]].sub(spdStates[1])
                # # sp1
                # print(temparr[b[0]].vel)
                # # sp2
                # print(temparr[b[1]].vel)

        # print(dataList.connectivity)

        # for a in range(len(dataList.data)):
        #     print(a)
        #     temparr.append(self.ambientSpace.acceleration(dataList.data[a].clone()).sub(spdStates[a]))
        return DataList(temparr, connectivity=dataList.connectivity)

    # //step forwards one timestep
    def step(self, dataList):

        # k1,k2,k3,k4;
        # temp;

        # //get the derivative
        k1 = self.deriveFunc(dataList)
        # print(k1.data[1].vel)
        temp1 = dataList.clone().add(k1.clone().multiplyScalar(self.stepSize/2.))
        # k1.multiplyScalar(self.ep)

        # //get k2
        # temp = dataList.clone().add(k1.clone().multiplyScalar(0.5))
        k2 = self.deriveFunc(temp1)
        temp2 = dataList.clone().add(k2.clone().multiplyScalar(self.stepSize/2.))
        # k2.multiplyScalar(self.stepSize)

        # //get k3
        # temp=dataList.clone().add(k2.clone().multiplyScalar(0.5))
        k3=self.deriveFunc(temp2)
        temp3 = dataList.clone().add(k2.clone().multiplyScalar(self.stepSize))
        # k3.multiplyScalar(self.stepSize)

        # //get k4
        # temp=dataList.clone().add(k3.multiplyScalar(1.))
        k4=self.deriveFunc(temp3)
        # k4.multiplyScalar(self.ep)

        # //add up results:
        total = k1 #//scale factor 1
        total.add(k2.multiplyScalar(2))
        total.add(k3.multiplyScalar(2))
        total.add(k4) #//scale factor 1
        total.multiplyScalar(self.stepSize/6.)

        # print(dataList.clone().updateBy(total).data[0].pos)
        # print(total.data[1].vel)
        # //move ahead one step
        return dataList.clone().updateBy(total)