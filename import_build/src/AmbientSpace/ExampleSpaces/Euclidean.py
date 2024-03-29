import numpy as np

from src.AmbientSpace.Components.Geometry import Geometry
from src.AmbientSpace.Components.Model import Model
from src.AmbientSpace.Components.Obstacle import Obstacle

from src.AmbientSpace.AmbientSpace import AmbientSpace

from src.Computation.State import State

# // -------------------------------------------------------------
# // 3-Dimensional Euclidean Space Information
# // -------------------------------------------------------------

# // -------------------------------------------------------------
# // Geometry Information
# // -------------------------------------------------------------

def eucMetricTensor(pos):
    return np.identity(3)

def eucChristoffel(state):
    return np.zeros(3)

def eucDistance(pos1, pos2):
    return np.sqrt(np.dot(np.subtract(pos1.copy(),pos2.copy()),np.subtract(pos1.copy(),pos2.copy())))


eucSpace = Geometry(
    eucMetricTensor,
    eucChristoffel,
    eucDistance
    )

# // -------------------------------------------------------------
# // Model Information
# // -------------------------------------------------------------

def identityR3(coords):
    return coords

def unitScaling(pos):
    return 1.

eucModel = Model(identityR3,unitScaling)

# // -------------------------------------------------------------
# // Obstacle/Bounding Ball Information
# // -------------------------------------------------------------

# // Sphere Bounding Box

# Default bounding box size (radius)
R = 6.

def distToSphere(pos):
    return R-np.sqrt(pos[0]**2. + pos[1]**2. + pos[2]**2.)

sphereObstacle = Obstacle(
    distToSphere,
    R
)

# // Cube Bounding Box

# Default bounding box size (cube side length)
boxSize = 6.

def distToBox(pos):
    xWall = boxSize - abs(pos.x)
    yWall = boxSize - abs(pos.y)
    zWall = boxSize - abs(pos.z)

    return min(xWall,min(yWall,zWall))

boxObstacle = Obstacle(
    distToBox,
    boxSize
)

# //package stuff up for export
euclidean = AmbientSpace( eucSpace, eucModel, sphereObstacle)

