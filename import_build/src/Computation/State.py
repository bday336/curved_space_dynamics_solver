# //The class "State" stores a tangent vector to configuration space
# //in practice, this will the configuration space of a single ball system,
# //meaning that state.pos will be the ball's center and state.vel will be
# //the balls velocity.  But the attributes can be anything,
# //so long as they implement the following methods:


# // .pos and .vel need to implement
# // .clone()
# // .vel needs to implement vector space operations:
# // .add(), .sub() . multiplyScalar(), .applyMatrix3()


class State:

    # //build a state from the input of an object storing position data
    # //and an object storing velocity data
    def __init__(self, pos, vel):
        self.pos = pos.copy()
        self.vel = vel.copy()

    # //make a copy of a given state (not just reference it in memory)
    # //and return the copy
    def clone(self):
        return State(self.pos.copy(), self.vel.copy())


    # //add the velocity of a given state to the current
    def add(self, state ):
        self.vel = self.vel + state.vel

    # //subtract the velocity of a given state from the current
    def sub(self, state ):
        self.vel = self.vel - state.vel

    # //scale the velocity of the current state by a factor
    def multiplyScalar(self, k ):
        self.vel = self.vel * k

    # //take the directional derivative of a function fn
    # // at pos in direction vel:
    def differentiate(self,fn):

        eps = 0.00001
        pos1 = self.pos.copy().add(self.vel.copy().multiplyScalar(-eps/2))
        pos2 = self.pos.copy().add(self.vel.copy().multiplyScalar(eps/2))

        dval = fn(pos2)-fn(pos1)
        return  dval/eps



    # //move a state infintesimally along its tangent direction
    def flow(self, ep):
        self.pos = self.pos + self.vel.copy() * ep

    # //update a state (a tangent vector) by infinitesimally flowing along a
    # //differential to the state: a pair dState of a velocity and acceleration
    def updateBy(self, dState ):
        self.pos = self.pos + dState.vel
        self.vel = self.vel + dState.acc