import numpy as np

import src.Computation.State as State

# A class used to store information about obstacle in space (e.i. bounding box)
class Obstacle:
    """
    A class used to store information about obstacle in space (e.i. bounding box)

    ...

    Attributes
    ----------
    distance : function
        function to calculate the distance between two points in the ambient space

    size : float
        Rough size of the bouding box of obstacle in geometric units

    """

    def __init__(self, distance, size):
        self.distance=distance
        self.size=size

