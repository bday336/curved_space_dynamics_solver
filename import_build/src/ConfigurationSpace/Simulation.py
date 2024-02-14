
# import {DataList} from "../Computation/DataList.js";
# import {RungeKutta} from "../Computation/RungeKutta.js";



# //RIDICULOUS: need to figure out how not to have to import these :(
# import{ ambientSpace, configurationSpace } from "../setup.js";

from src.Computation.DataList import DataList
from src.Computation.RungeKutta import RungeKutta

class Simulation:

    def __init__(self, ambientSpace, configurationSpace, dataList, sysparam, stepSize):
        self.ambientSpace = ambientSpace
        self.configurationSpace = configurationSpace

        self.dataList = dataList
        self.sysparam = sysparam
        self.stepSize = stepSize

        # //to set when intersecting
        self.ball_collisions = [] 
        self.obstacle_collisions = []

        self.data_container = []
        self.data_container.append(self.dataList)


        # //build an integrator
        # //get the function which takes the derivative of each element of a stateList:
        # //using ambientSpace.acceleration will allow us to use external potentials without changing code
        # self.derive = self.deriveFunc()

        self.integrator = RungeKutta(self.ambientSpace,self.stepSize)

    def deriveFunc(self, dataList):
        temparr = []
        for a in dataList.data:
            temparr.append(self.ambientSpace.acceleration(a.clone()))
        return DataList(temparr)


    def detectCollision(self):

        self.obstacle_collisions = self.configurationSpace.obstacleCollisions(self.dataList)
        self.ball_collisions = self.configurationSpace.ballCollisions(self.dataList)

        # print("obstacle detection")
        # print(self.obstacle_collisions)

        # print("ball detection")
        # print(self.ball_collisions)

        # Collision detected
        if len(self.obstacle_collisions) != 0 or len(self.ball_collisions) != 0:
            return True
        # No collision detected
        else:
            return False


    def smoothDynamics(self):
        # print(self.dataList)
        self.dataList = self.integrator.step(self.dataList)
        # print(self.dataList)

    def collisionDynamics(self):

        # //get the normal vector to the boundary of configuration space
        bdyNormal = self.configurationSpace.boundaryNormal(
            self.dataList,
            self.ball_collisions,
            self.obstacle_collisions
        )

        # //update the state by reflecting in the boundary normal
        # //with respect to the configuration space's metric
        self.dataList = self.configurationSpace.reflectIn(self.dataList, bdyNormal)

    def step(self):

        # //get the points of collision, if there are any
        collide = self.detectCollision()
        # print(collide)

        if collide :
            self.collisionDynamics()
            self.data_container.append(self.dataList)

        else:
            # //then after they've been resolved, run smooth dynamics
            self.smoothDynamics()
            self.data_container.append(self.dataList)

