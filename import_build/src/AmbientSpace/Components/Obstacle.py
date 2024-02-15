import numpy as np

import src.Computation.State as State



class Obstacle:
    def __init__(self, distance, geometry, size=None, generateState=None):
        self.distance=distance;
        self.geometry=geometry;


        # //size gives the rough size of the bounding box IN GEOMETRIC UNITS: useful for setting initial conditions
        # //(ie choosing radii)
        if(size!=None):
            self.size=size
        else:
            self.size = 1.


        # //generateState makes a random state IN COORDINATES that is not hitting the obstacle
        # if(generateState != None):
        #     self.generateState = generateState
        # else:
        #     self.generateState =  State(randomVec3Ball(1),randomVec3Ball(1))

