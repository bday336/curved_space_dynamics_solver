from src.Computation.DataList import DataList
from src.Computation.RungeKutta import RungeKutta
from src.Computation.Gauss1 import Gauss1
from src.Computation.Gauss1test import Gauss1test
from src.Computation.Gauss2 import Gauss2
from src.Computation.Gauss2test import Gauss2test
from src.Computation.Gauss3 import Gauss3
from src.Computation.Radau2 import Radau2
from src.Computation.Radau2test import Radau2test
from src.Computation.RigidRadau2 import RigidRadau2
from src.Computation.Radau3 import Radau3
from src.Computation.RigidRadau3 import RigidRadau3

class Simulation:
    """
    A class used to store information about simulation

    ...

    Attributes
    ----------
    ambientSpace : object
        AmbientSpace object describing ambient space containing system

    dataList : object
        DataList object describing system (collection of vertex information)

    configurationSpace : object
        ConfigurationSpace object describing configuration space of simulation systems

    stepSize : float
        time step for each integration step

    Methods
    -------
    detectCollision()
        Check to determine if any collisions have occurred between vertices or between vertices and obstacles
        Returns boolean

    smoothDynamics()
        Perform one step of numerical integration
        Updates self.dataList to post integration values

    collisionDynamics()
        Calculate boundary normal in configuration space and update self.dataList to post collision values
        Updates self.dataList to post collision values

    step()
        Perform one step of simulation:
        1) Detect if collision present
        2) if collision present run collsionDyamics()
        3) if no collision present run smoothDynamics()
    """

    def __init__(self, ambientSpace, dataList, configurationSpace, stepSize, solver_method):
        self.ambientSpace = ambientSpace
        self.dataList = dataList
        self.configurationSpace = configurationSpace

        self.solver_method = solver_method
        self.stepSize = stepSize

        # //to set when intersecting
        self.ball_collisions = [] 
        self.obstacle_collisions = []

        # Container for all simulation data
        self.data_container = []
        self.data_container.append(self.dataList)


        # //build an integrator
        # //get the function which takes the derivative of each element of a stateList:
        # //using ambientSpace.acceleration will allow us to use external potentials without changing code
        # self.derive = self.deriveFunc()

        if self.solver_method == "RungeKutta":
            self.integrator = RungeKutta(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Gauss1":
            self.integrator = Gauss1(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Gauss1test":
            self.integrator = Gauss1test(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Gauss2":
            self.integrator = Gauss2(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Gauss2test":
            self.integrator = Gauss2test(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Gauss3":
            self.integrator = Gauss3(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Radau2":
            self.integrator = Radau2(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Radau2test":
            self.integrator = Radau2test(self.ambientSpace,self.stepSize)
        elif self.solver_method == "RigidRadau2":
            self.integrator = RigidRadau2(self.ambientSpace,self.stepSize)
        elif self.solver_method == "Radau3":
            self.integrator = Radau3(self.ambientSpace,self.stepSize)
        elif self.solver_method == "RigidRadau3":
            self.integrator = RigidRadau3(self.ambientSpace,self.stepSize)

    # def deriveFunc(self, dataList):
    #     temparr = []
    #     for a in dataList.data:
    #         temparr.append(self.ambientSpace.acceleration(a.clone()))
    #     return DataList(temparr)


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

