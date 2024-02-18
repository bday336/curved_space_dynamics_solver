# //The class "dState" stores a differential of a
# // tangent vector to configuration space
# //in practice, this will the configuration space of a single ball system,
# //meaning that state.vel will be the ball's velocity and state.acc will be
# //covariant derivative with respect to the metric on the ambient space
# //so long as they implement the following methods:


# // .vel and .acc need to implement
# // .clone(), as well as the  vector space operations:
# // .add(), .sub() . multiplyScalar(), .applyMatrix3()

import numpy as np


class dState:
    """
    A class used to store information about dstate
    Describes vector

    ...

    Attributes
    ----------
    vel : array
        array of velocity data

    acc : array
        array of acceleration data

    Methods
    -------
    clone()
        Generate copy of self
        Returns dState clone

    add(dstate)
        Add velocity and acceleration of dstate to self.vel and self.acc, respectively
        Returns self

    sub(dstate)
        Subtract velocity and acceleration of dstate from self.vel and self.acc, respectively
        Returns self

    multiplyScalar(k)
        Scale self.vel and self.acc by scalar value k
        Returns self
    """

    # //build a dState from the input of an object storing velocity data
    # //and an object storing acceleration data
    def __init__(self, vel, acc):
        self.vel=vel.copy()
        self.acc=acc.copy()

    # //make a copy of a given dState (not just reference it in memory)
    # //and return the copy
    def clone(self):
        return  dState(self.vel.copy(),self.acc.copy())


    # //add the velocity AND of a given state to the current
    def add(self, dState):
        self.vel = np.add(self.vel, dState.vel)
        self.acc = np.add(self.acc, dState.acc)
        return self

    # //subtract the velocity AND acceleration of a given state from the current
    def sub(self, dState ):
        self.vel = np.subtract(self.vel, dState.vel )
        self.acc = np.subtract(self.acc, dState.acc )
        return self


    # //scale the velocity AND acceleration of the current state by a factor
    def multiplyScalar (self, k ):
        self.vel = self.vel * k
        self.acc = self.acc * k 
        return self
