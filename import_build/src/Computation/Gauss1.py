# // a class for integrating equations of motion called "integrator" and one specific implementation, Runge Kutta
# //derive is a function taking a state to state (now storing velocity and acceleration instead of position and velocity)
# //items fed into RungeKutta need to have the following methods available:
# //.add, .multiplyScalar, .clone

import numpy as np
from numpy import sin, cos, tan, cosh, sinh, tanh, arccosh, sqrt, array, zeros

from src.Computation.DataList import DataList
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

    # Jacobian will be constructed to be largely diagonal to help with efficiently for larger (sparse matrix) systems - block diagonal jacobian
    # This is in contrast with previous iteration of the solver design where top half was velocity residual expressions and bottom half were acceleration residual expressions
    def gen_sys_jacobian(self, dataList):
        # Dimension of ambient space (hard coded to 3D at the moment)
        ambient_dim = 3
        # Number of vertices in system
        vert_num = len(dataList.data)
        # Initialize jacobian matrix (multiply by 2 since considering residuals of velocity and accelerations for each vertex)
        jacobian = np.zeros((2*ambient_dim*vert_num,2*ambient_dim*vert_num))

        # Populate vertex data (block diagonal information in jacobian)
        for a in range(vert_num):
            tempmat = self.ambientSpace.geometry.funcDict["vertex_jac"](dataList.data[a])
            jacobian[a*2*ambient_dim : a*2*ambient_dim + 2 * ambient_dim * vert_num , a*2*ambient_dim : a*2*ambient_dim + 2 * ambient_dim * vert_num] = tempmat.copy()

        ks,x = [1,1]
        m = 1.

    # a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    # v,ks,x,m = params

    # # Distance function
    # d12 = D12(a1, b1, g1, a2, b2, g2)
    # # First derivatives of distance function
    # da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    # db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    # dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    # da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    # db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    # dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    # # Second derivatives of distance function
    # da1d12a1=da1D12a1(a1, b1, g1, a2, b2, g2)
    # db1d12a1=db1D12a1(a1, b1, g1, a2, b2, g2)
    # dg1d12a1=dg1D12a1(a1, b1, g1, a2, b2, g2)
    # da2d12a1=da2D12a1(a1, b1, g1, a2, b2, g2)
    # db2d12a1=db2D12a1(a1, b1, g1, a2, b2, g2)
    # dg2d12a1=dg2D12a1(a1, b1, g1, a2, b2, g2)
    
    # da1d12b1=db1d12a1
    # db1d12b1=db1D12b1(a1, b1, g1, a2, b2, g2)
    # dg1d12b1=dg1D12b1(a1, b1, g1, a2, b2, g2)
    # da2d12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)
    # db2d12b1=db2D12b1(a1, b1, g1, a2, b2, g2)
    # dg2d12b1=dg2D12b1(a1, b1, g1, a2, b2, g2)

    # da1d12g1=dg1d12a1
    # db1d12g1=dg1d12b1
    # dg1d12g1=dg1D12g1(a1, b1, g1, a2, b2, g2)
    # da2d12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
    # db2d12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)
    # dg2d12g1=dg2D12g1(a1, b1, g1, a2, b2, g2)

    # da1d12a2=da2d12a1
    # db1d12a2=da2d12b1
    # dg1d12a2=da2d12g1
    # da2d12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
    # db2d12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
    # dg2d12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

    # da1d12b2=db2d12a1
    # db1d12b2=db2d12b1
    # dg1d12b2=db2d12g1
    # da2d12b2=db2d12a2
    # db2d12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
    # dg2d12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

    # da1d12g2=dg2d12a1
    # db1d12g2=dg2d12b1
    # dg1d12g2=dg2d12g1
    # da2d12g2=dg2d12a2
    # db2d12g2=dg2d12b2
    # dg2d12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)

    # #---------- particle 1

    # da1riga1 = cosh(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + da2da1V12(m, 1., ks, x, d12, da1d12, da1d12, 0., da1d12a1)
    # db1riga1 = cos(b1)*sin(b1)*sinh(2.*a1)*gd1**2.           + da2da1V12(m, 1., ks, x, d12, da1d12, db1d12, 0., db1d12a1)
    # dg1riga1 = 0.                                            + da2da1V12(m, 1., ks, x, d12, da1d12, dg1d12, 0., dg1d12a1)

    # da2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, da2d12, 0., da2d12a1)
    # db2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, db2d12, 0., db2d12a1)
    # dg2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, dg2d12, 0., dg2d12a1)

    # dad1riga1 = 0.                          + 0.
    # dbd1riga1 = sinh(2.*a1)*bd1             + 0.
    # dgd1riga1 = sin(b1)**2.*sinh(2.*a1)*gd1 + 0.

    # dad2riga1 = 0. + 0.
    # dbd2riga1 = 0. + 0.
    # dgd2riga1 = 0. + 0.

    # #----------

    # da1rigb1 = 2./sinh(a1)**2.*(ad1*bd1) + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, da1d12, sinh(2.*a1), da1d12b1)
    # db1rigb1 = cos(2.*b1)*gd1**2.        + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, db1d12, 0.,          db1d12b1)
    # dg1rigb1 = 0.                        + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, dg1d12, 0.,          dg1d12b1)

    # da2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, da2d12, 0., da2d12b1)
    # db2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, db2d12, 0., db2d12b1)
    # dg2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, dg2d12, 0., dg2d12b1)

    # dad1rigb1 = -2./tanh(a1)*bd1 + 0.
    # dbd1rigb1 = -2./tanh(a1)*ad1 + 0.
    # dgd1rigb1 = sin(2.*b1)*g1    + 0.

    # dad2rigb1 = 0. + 0.
    # dbd2rigb1 = 0. + 0.
    # dgd2rigb1 = 0. + 0.

    # #------------

    # da1rigg1 = 2./sinh(a1)**2.*(ad1*gd1) + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da1d12, sinh(2.*a1)*sin(b1)**2., da1d12g1)
    # db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)  + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db1d12, sin(2.*b1)*sinh(a1)**2., db1d12g1)
    # dg1rigg1 = 0.                        + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg1d12, 0.,                      dg1d12g1)

    # da2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da2d12, 0., da2d12g1)
    # db2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db2d12, 0., db2d12g1)
    # dg2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg2d12, 0., dg2d12g1)

    # dad1rigg1 = -2./tanh(a1)*gd1                 + 0.
    # dbd1rigg1 = -2./tan(b1)*gd1                  + 0.
    # dgd1rigg1 = -2.*(ad1/tanh(a1) + bd1/tan(b1)) + 0.

    # dad2rigg1 = 0. + 0.
    # dbd2rigg1 = 0. + 0.
    # dgd2rigg1 = 0. + 0.

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

                # Incorporate into correct element of jacobian matrix

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

    def dynfunc_h3simgeo(self,state_vec,params):
        
        a1,b1,g1,ad1,bd1,gd1 = state_vec
        v,m = params


        riga1 = 1./2.*sinh(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) + 0.
        rigb1 = -2. * ad1*bd1/tanh(a1) + .5*sin(2.*b1)*gd1**2         + 0.
        rigg1 = -2. * ad1*gd1/tanh(a1) - 2.*bd1*gd1/tan(b1)           + 0.

        return np.array(
            [ad1,bd1,gd1,
            riga1,rigb1,rigg1])

    def difffuncgauss1s(self, startvec, params, dynfunc, k1, dt):
        # Set to run Gauss 1-stage method
        a11 = 1./2.

        return array([
            k1 - dynfunc(startvec + (a11*k1)*dt, params)
        ]).flatten()

    def gausss1(self, startvec, params, dynfunc, dynjac, dt, tol = 1e-15, imax = 100):
        # Set to run Gauss 1-stage method
        a11 = 1./2.
        bs1 = 1.

        # Initial Guess - Explicit Euler
        k = dynfunc(startvec, params)
        x1guess = startvec + (1./2.)*dt*k
        k1 = dynfunc(x1guess, params)

        # Check Error Before iterations
        er = self.difffuncgauss1s(startvec, params, dynfunc, k1, dt)

        # Begin Iterations
        counter = 0
        while (np.linalg.norm(er) >= tol and counter <= imax):
            j1 = dynjac(startvec + (a11*k1)*dt, params)
            
            fullj = np.block([
                [np.eye(k.shape[0]) - dt*a11*j1]
            ])

            linsolve = np.linalg.solve(fullj,-er)

            k1 = k1 + linsolve[0:k.shape[0]]

            er = self.difffuncgauss1s(startvec, params, dynfunc, k1, dt)

            counter += 1

        startvec = startvec + dt*(bs1*k1)
        return startvec.copy()

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