# // a class for integrating equations of motion called "integrator" and one specific implementation, Runge Kutta
# //derive is a function taking a state to state (now storing velocity and acceleration instead of position and velocity)
# //items fed into RungeKutta need to have the following methods available:
# //.add, .multiplyScalar, .clone

from src.Computation.DataList import DataList

# //implementing the Rk4 Scheme for arbitrary classes that have clone add and multiplyScalar
# //will use this possibly on individual states, or on entire DataLists!
class RungeKutta:

    def __init__(self, ambientSpace, stepSize):
        # self.derive=derive
        self.ambientSpace = ambientSpace
        self.stepSize = stepSize

    def deriveFunc(self, dataList):
        temparr = []
        for a in dataList.data:
            temparr.append(self.ambientSpace.acceleration(a.clone()))
        return DataList(temparr)

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