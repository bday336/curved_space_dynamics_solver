import numpy as np

from src.AmbientSpace.Components.Geometry import Geometry
from src.AmbientSpace.Components.Model import Model
from src.AmbientSpace.Components.Obstacle import Obstacle

from src.AmbientSpace.AmbientSpace import AmbientSpace

from src.Computation.State import State

# // -------------------------------------------------------------
# // Example Space Information Template
# // -------------------------------------------------------------

# // -------------------------------------------------------------
# // Geometry Information
# // -------------------------------------------------------------

def MetricTensor(pos):
    return 

def Christoffel(state):
    return 

def Distance(pos1, pos2):
    return 


space = Geometry(
    MetricTensor,
    Christoffel,
    Distance
    )

# // -------------------------------------------------------------
# // Model Information
# // -------------------------------------------------------------

def identityR3(coords):
    return 

def unitScaling(pos):
    return 

model = Model(identityR3,unitScaling)

# // -------------------------------------------------------------
# // Obstacle/Bounding Ball Information
# // -------------------------------------------------------------

# // Sphere Bounding Box

# Default bounding box size (radius)
R = 6.

def distToSphere(pos):
    return 

sphereObstacle = Obstacle(
    distToSphere,
    R
)

# //package stuff up for export
templateSpace = AmbientSpace( space, model, sphereObstacle)

