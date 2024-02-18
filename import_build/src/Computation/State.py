# //The class "State" stores a tangent vector to configuration space
# //in practice, this will the configuration space of a single ball system,
# //meaning that state.pos will be the ball's center and state.vel will be
# //the balls velocity.  But the attributes can be anything,
# //so long as they implement the following methods:


# // .pos and .vel need to implement
# // .clone()
# // .vel needs to implement vector space operations:
# // .add(), .sub() . multiplyScalar(), .applyMatrix3()

import numpy as np

class State:
    """
    A class used to store information about state (vertices of system)
    Describes tangent vector at vertex (velocity of vertex)

    ...

    Attributes
    ----------
    pos : array
        array of position of state in the ambient space

    vel : array
        array of velocity of state in the ambient space

    Methods
    -------
    clone()
        Generate copy of self
        Returns State clone

    add(state)
        Add velocity of state to self.vel
        Returns self

    sub(state)
        Subtract velocity of state from self.vel
        Returns self

    multiplyScalar(k)
        Scale self.vel by scalar value k
        Returns self

    differentiate(fn)
        Calculate directional derivative of function fn at position self.pos in direction self.vel
        Returns scaled output of fn

    flow(eps)
        Move state infinitesimally (i.e. by eps) along its tangent vector
        Returns self

    updateBy(dState)
        Update state by infinitesimally flowing along a differential to the state (dState)
        Returns self
    """

    # //build a state from the input of an object storing position data
    # //and an object storing velocity data
    def __init__(self, pos, vel):
        self.pos=pos.copy()
        self.vel=vel.copy()

    # //make a copy of a given state (not just reference it in memory)
    # //and return the copy
    def clone(self):
        return  State(self.pos.copy(), self.vel.copy())


    # //add the velocity of a given state to the current
    def add(self, state ):
        # print(self.vel)
        self.vel = np.add(self.vel,state.vel)
        return self

    # //subtract the velocity of a given state from the current
    def sub(self, state ):
        self.vel = np.subtract(self.vel,state.vel)
        return self

    # //scale the velocity of the current state by a factor
    def multiplyScalar(self, k ):
        self.vel = self.vel * k
        return self

    # //take the directional derivative of a function fn
    # // at pos in direction vel:
    def differentiate(self,fn):

        eps = 0.00001
        pos1 = np.add(self.pos.copy(),self.vel.copy()*-eps/2)
        pos2 = np.add(self.pos.copy(),self.vel.copy()*eps/2)

        dval = fn(pos2)-fn(pos1)
        return  dval/eps

    # //move a state infintesimally along its tangent direction
    def flow(self, eps):
        self.pos = np.add(self.pos,self.vel.copy()*eps)
        return self

    # //update a state (a tangent vector) by infinitesimally flowing along a
    # //differential to the state: a pair dState of a velocity and acceleration
    def updateBy(self, dState ):
        self.pos = np.add(self.pos,dState.vel)
        self.vel = np.add(self.vel,dState.acc)
        return self