import numpy as np

from src.AmbientSpace.Components.Geometry import Geometry
from src.AmbientSpace.Components.Model import Model
from src.AmbientSpace.Components.Obstacle import Obstacle

from src.AmbientSpace.AmbientSpace import AmbientSpace

from src.Computation.State import State

# // -------------------------------------------------------------
# // 3-Dimensional Hyperbolic Space Information (Poincare Ball Model)
# // -------------------------------------------------------------

# // -------------------------------------------------------------
# // Helper Functions
# // -------------------------------------------------------------

def minkowskiDot(u,v):
    return u[3]*v[3] - ( u[0]*v[0] + u[1]*v[1] + u[2]*v[2] )

# //distance on hyperboloid:
def hyperboloidDistance(u,v):
    return np.arccosh(abs(minkowskiDot(u,v)))

# //map from poincare ball to the hyperboloid:
def toHyperboloid(pos):

    len2 = pos[0]**2. + pos[1]**2. + pos[2]**2.
    w = 1 + len2
    p = np.array([2.*pos[0],2.*pos[1],2.* pos[2],w]) / (1-len2)

    return p

# // -------------------------------------------------------------
# // Geometry Information
# // -------------------------------------------------------------


# //this metric is conformal to the Euclidean plane
# //so, its defined by a conformal factor
def conformalFactor(pos):

    r2 = pos[0]**2. + pos[1]**2. + pos[2]**2.
    diff = 1-r2
    diff2 = diff*diff

    return  4./(diff2)

def hypMetricTensor(pos):

    # //just multiply the identity by the conformal factor
    scale = conformalFactor(pos)
    return scale * np.identity(3)

def hypChristoffel(state):

    pos = state.pos.copy()
    x = pos[0]
    y = pos[1]
    z = pos[2]

    vel = state.vel.copy()
    xP = vel[0]
    yP = vel[1]
    zP = vel[2]

    xP2 = xP*xP
    yP2 = yP*yP
    zP2 = zP*zP

    denom = pos[0]**2. + pos[1]**2. + pos[2]**2. - 1.

    xPP = 2*x * ( xP2 - yP2 - zP2 ) + 4 * xP * ( y*yP + z*zP )
    yPP = 2*y * ( yP2 - xP2 - zP2 ) + 4 * yP * ( x*xP + z*zP )
    zPP = 2*z * ( zP2 - xP2 - yP2 ) + 4 * zP * ( x*xP + y*yP )

    acc =  np.array([xPP,yPP,zPP]) / denom

    return acc


def hypDistance(pos1,pos2):

    u = toHyperboloid(pos1)
    v = toHyperboloid(pos2)
    return hyperboloidDistance(u,v)


hypSpace = Geometry(
    hypMetricTensor,
    hypChristoffel,
    hypDistance
    )


# // -------------------------------------------------------------
# // Model Information
# // -------------------------------------------------------------

# //there is no model for this space: its a metric directly on R3!
# //though, its all drawn inside of a unit ball: so let's scale it up

zoom = 6.

def hypCoordsToModel(coords):
    return coords.copy().multiplyScalar(zoom)


# //the scaling factor is computed from the metric tensor:
# //this metric tensor is conformal so its easy: sqrt(conformalCoef)

def hypModelScaling(modelPos):

    # //unscale position back to true Poincare ball:
    coordPos = modelPos.copy() / zoom
    scale = conformalFactor(coordPos)
    return zoom/np.sqrt(scale)


hypModel = Model(hypCoordsToModel, hypModelScaling)




# // -------------------------------------------------------------
# // Obstacle/Bounding Ball Information
# // -------------------------------------------------------------

# //a sphere of radius R
coordSize = 0.8
sphereSize = hypDistance(np.array([0,0,0]), np.array([coordSize,0,0]))

# //the metric distance to the origin of the poincare disk is the arctanh of the norm:
def distToSphere(pos):
    center = np.array([0,0,0])
    dist =  hypDistance(pos,center)
    return sphereSize - dist

sphereObstacle = Obstacle(
    distToSphere,
    sphereSize
)

# //package stuff up for export
hyperbolic = AmbientSpace( hypSpace, hypModel, sphereObstacle)
