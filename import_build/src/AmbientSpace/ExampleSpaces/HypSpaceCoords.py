# // import {
# //     Matrix3,
# //     SphereBufferGeometry,
# //     BoxBufferGeometry,
# //     Vector3,
# //     Vector4
# // } from "../../../3party/three/build/three.module.js";

# // import { Geometry } from "../Components/Geometry.js";
# // import {Model} from "../Components/Model.js";
# // import {Obstacle} from "../Components/Obstacle.js";

# // import {AmbientSpace} from "../AmbientSpace.js";

import numpy as np

from src.AmbientSpace.Components.Geometry import Geometry
from src.AmbientSpace.Components.Model import Model
from src.AmbientSpace.Components.Obstacle import Obstacle

from src.AmbientSpace.AmbientSpace import AmbientSpace

from src.Computation.State import State

# // -------------------------------------------------------------
# //some geometry stuff:
# // -------------------------------------------------------------


def minkowskiDot(u,v):
    return u[3]*v[3] - ( u[0]*v[0] + u[1]*v[1] + u[2]*v[2] )


# //distance on hyperboloid:
def hyperboloidDistance(u,v):
    # print("check dotproduct")
    # print(abs(minkowskiDot(u,v)))
    return np.arccosh(abs(minkowskiDot(u,v)))


# //coordinates mapping onto the hyperboloid model of H3
def coords(pos):

    alpha = pos[0]
    beta  = pos[1]
    gamma = pos[2]

# Translational parameterization
    # x = np.sinh(alpha)
    # y = np.cosh(alpha)*np.sinh(beta)
    # z = np.cosh(alpha)*np.cosh(beta)*np.sinh(gamma)
    # w = np.cosh(alpha)*np.cosh(beta)*np.cosh(gamma)

# Rotational parameterization
    x = np.sinh(alpha)*np.sin(beta)*np.cos(gamma)
    y = np.sinh(alpha)*np.sin(beta)*np.sin(gamma)
    z = np.sinh(alpha)*np.cos(beta)
    w = np.cosh(alpha)

    return np.array([x,y,z,w])



def hypMetricTensor(pos):

    alpha = pos[0]
    beta  = pos[1]
    gamma = pos[2]

# Translational parameterization
    # cosh2alpha = np.cosh(alpha)*np.cosh(alpha)
    # cosh2beta  = np.cosh(beta)*np.cosh(beta)

    # g11 = 1.
    # g22 = cosh2alpha
    # g33 = cosh2alpha*cosh2beta

# Rotational parameterization
    sinh2alpha = np.sinh(alpha) * np.sinh(alpha)
    sinh2beta  = np.sinh(beta) *  np.sinh(beta)

    g11 = 1.
    g22 = sinh2alpha
    g33 = sinh2alpha*sinh2beta

    return np.array([
        [g11,0,0],
        [0,g22,0],
        [0,0,g33]
    ])


def hypDistance(pos1, pos2):

    u = coords(pos1)
    v = coords(pos2)

    return hyperboloidDistance(u,v)

def hypChristoffel(state):

    pos = state.pos.copy()
    vel = state.vel.copy()

    alpha = pos[0]
    beta  = pos[1]
    gamma = pos[2]

    dalpha = vel[0]
    dbeta  = vel[1]
    dgamma = vel[2]
    
# Translational Parameterization
    # ddalpha = np.cosh(alpha)*np.sinh(alpha)*(dbeta*dbeta+np.cosh(beta)*np.cosh(beta)*dgamma*dgamma)
    # ddbeta  = np.cosh(beta)*np.sinh(beta)*dgamma*dgamma - 2*np.tanh(alpha)*dbeta*dalpha
    # ddgamma = -2*dgamma*(dbeta*np.tanh(beta)+dalpha*np.tanh(alpha))

# Rotational Parameterization
    ddalpha = .5 * np.sinh(2. * alpha) * (dbeta**2. + np.sin(beta)**2. * dgamma**2.)
    ddbeta  = np.sin(2. * beta) * dgamma**2. - 2. / np.tanh(alpha) * dbeta * dalpha
    ddgamma = -2. * dgamma * (dbeta / np.tan(beta) + dalpha / np.tanh(alpha))

    acc = np.array([ddalpha, ddbeta, ddgamma])

    return acc


hypSpace = Geometry(
    hypMetricTensor,
    hypChristoffel,
    hypDistance
    )





# // -------------------------------------------------------------
# //model of Euclidean space : do nothing
# // -------------------------------------------------------------

# //scale of the model on the screen
modelScale = 6.


def toPoincareBall(coord):

    p = coords(coord);

    x, y, z, w = p

    scale = modelScale / (1. + w)

    return np.array([x,y,z]) * scale



def pbScaling(pos):
    len2 = pos[0]**2. + pos[1]**2. + pos[2]**2.
    scale = modelScale**2. - len2
    return 4. * scale / (modelScale**2.)


hypModel = Model(toPoincareBall,pbScaling)




# // -------------------------------------------------------------
# //obstacles to contain balls in Euclidean Space
# // -------------------------------------------------------------

# //a sphere of a fixed radius
obstacleSize =2.

def distToSphere(pos):
    # //center point of H3 in coordinates:
    center = np.zeros(3)
    # //distance from center to position
    dist = hypDistance(pos.copy(),center)
    # //how far is this from the boundary sphere of radius 5?
    return obstacleSize - dist


rad = modelScale * np.tanh(obstacleSize)

sphereGeom = None #new SphereBufferGeometry(rad, 64, 32);

sphereObstacle = Obstacle(
    distToSphere,
    sphereGeom,
    obstacleSize
)







# //package stuff up for export
hyperbolic = AmbientSpace( hypSpace, hypModel, sphereObstacle)

# export { hyperbolic };


