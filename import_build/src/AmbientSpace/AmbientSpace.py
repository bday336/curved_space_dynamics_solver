# The class representing the ambient space of the geometry
# implementing what we have in the paper (except; I've left out potential energy for now, easy to come back and add)
# will just require a .potential()
# and an update to .acceleration() to not just be geodesic acceleration but also use this.gradient(this.potential(pos))

class AmbientSpace:
    def __init__(self, geometry, obstacle):
        self.geometry = geometry
        self.obstacle = obstacle

    # For adding force field
    def acceleration(self, state):
        return self.geometry.covariantAcceleration(state)

    # Dot product of geometry
    def dot(self, state1, state2):
        return self.geometry.dot(state1,state2)

    # Gradient with respect to geometry
    def gradient(self, fn, pos):
        return self.geometry.gradient(fn,pos)

    # Distance function of geometry
    def distance(self, pos1, pos2):
        return self.geometry.distance(pos1, pos2)

    # Distance function to the obstacle / bounding box of the simulation
    def distToObstacle(self, pos):
        return self.obstacle.distance(pos)