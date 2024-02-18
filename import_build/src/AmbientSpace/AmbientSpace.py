# The class representing the ambient space of the geometry
# implementing what we have in the paper (except; I've left out potential energy for now, easy to come back and add)
# will just require a .potential()
# and an update to .acceleration() to not just be geodesic acceleration but also use this.gradient(this.potential(pos))

class AmbientSpace:
    """
    A class used to store information about ambient space

    ...

    Attributes
    ----------
    geometry : object
        Geometry object describing ambient space of simulation

    model : object
        Model object describing how to visualize the data

    obstacle : object
        Obstacle object describing obstracle in ambient space

    Methods
    -------
    acceleration(state)
        Add force field to state
        Returns dState from covariant derivative of ambient space

    dot(state1, state2)
        Calculates the dot product of velocity vectors for specificed states
        Returns scalar value of dot product

    gradient(fn, pos)
        Calculate the gradient of function fn at position pos
        Returns differential (gradient) of function at given point

    toR3(pos)
        Project data using specified model for visualization
        Returns dictonary of visualization information

    distance(pos1, pos2)
        Calculates distance between position pos1 and pos2 in ambient space
        Returns scalar distance

    distanceToObstacle(pos)
        Calculates distance between obstacle/bounding box of simulation environment and position pos
        Returns scalar distance
    """

    def __init__(self, geometry, model, obstacle):
        self.geometry = geometry
        self.model    = model
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

    # Projection from model for visualization
    def toR3(self, pos):
        posR3 = self.model.toR3(pos)
        return {pos: posR3, scaling: self.model.relativeScaling(posR3)}

    # Distance function of geometry
    def distance(self, pos1, pos2):
        return self.geometry.distance(pos1, pos2)

    # Distance function to the obstacle / bounding box of the simulation
    def distToObstacle(self, pos):
        return self.obstacle.distance(pos)