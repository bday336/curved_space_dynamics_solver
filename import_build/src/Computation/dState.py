# //The class "dState" stores a differential of a
# // tangent vector to configuration space
# //in practice, this will the configuration space of a single ball system,
# //meaning that state.vel will be the ball's velocity and state.acc will be
# //covariant derivative with respect to the metric on the ambient space
# //so long as they implement the following methods:


# // .vel and .acc need to implement
# // .clone(), as well as the  vector space operations:
# // .add(), .sub() . multiplyScalar(), .applyMatrix3()




class dState:

    # //build a dState from the input of an object storing velocity data
    # //and an object storing acceleration data
    def __init__(self, vel, acc):
        self.vel = vel.copy()
        self.acc = acc.copy()

    # //make a copy of a given dState (not just reference it in memory)
    # //and return the copy
    def clone(self):
        return  dState(self.vel.copy(),self.acc.copy())


    # //add the velocity AND of a given state to the current
    def add(self, dState):
        self.vel = self.vel + dState.vel
        self.acc = self.acc + dState.acc

    # //subtract the velocity AND acceleration of a given state from the current
    def sub(self, dState ):
        self.vel = self.vel - dState.vel
        self.acc = self.acc - dState.acc


    # //scale the velocity AND acceleration of the current state by a factor
    def multiplyScalar (self, k ):
        self.vel = self.vel * k 
        self.acc = self.acc * k 