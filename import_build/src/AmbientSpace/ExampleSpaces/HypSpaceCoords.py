import numpy as np
from numpy import cos, sin, tan, cosh, sinh, tanh, sqrt, arccosh

from src.AmbientSpace.Components.Geometry import Geometry
from src.AmbientSpace.Components.Model import Model
from src.AmbientSpace.Components.Obstacle import Obstacle

from src.AmbientSpace.AmbientSpace import AmbientSpace

from src.Computation.State import State

# // -------------------------------------------------------------
# // 3-Dimensional Hyperbolic Space Information (Hyperboloid Model)
# // -------------------------------------------------------------

# // -------------------------------------------------------------
# // Helper Functions
# // -------------------------------------------------------------

def minkowskiDot(u,v):
    return u[3]*v[3] - ( u[0]*v[0] + u[1]*v[1] + u[2]*v[2] )

# //distance on hyperboloid:
def hyperboloidDistance(u,v):
    return np.arccosh(abs(minkowskiDot(u,v)))

# //coordinates mapping onto the hyperboloid model of H3 in 4D minkowski space
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

# // -------------------------------------------------------------
# // Geometry Information
# // -------------------------------------------------------------

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
    ddbeta  = .5 * np.sin(2. * beta) * dgamma**2. - 2. / np.tanh(alpha) * dbeta * dalpha
    ddgamma = -2. * dgamma * (dbeta / np.tan(beta) + dalpha / np.tanh(alpha))

    acc = np.array([ddalpha, ddbeta, ddgamma])

    return acc

def hyp_gen_jacobian(state):

    pos = state.pos.copy()
    vel = state.vel.copy()

    alpha = pos[0]
    beta  = pos[1]
    gamma = pos[2]

    dalpha = vel[0]
    dbeta  = vel[1]
    dgamma = vel[2]

    # Below the labeling corresponds with:
    # daddalpha  -> derivative of ddalpha (see hypChristoffel()) with respect to alpha
    # ddaddalpha -> derivative of ddalpha (see hypChristoffel()) with respect to dalpha
    #---------- 

    daddalpha = cosh(2.*alpha) * (dbeta**2. + sin(beta)**2.*dgamma**2.)
    dbddalpha = cos(beta)*sin(beta)*sinh(2.*alpha)*dgamma**2.
    dgddalpha = 0.

    ddaddalpha = 0.
    ddbddalpha = sinh(2.*alpha)*dbeta
    ddgddalpha = sin(beta)**2.*sinh(2.*alpha)*dgamma

    #----------

    daddbeta = 2./sinh(alpha)**2.*(dalpha*dbeta)
    dbddbeta = cos(2.*beta)*dgamma**2.
    dgddbeta = 0.

    ddaddbeta = -2./tanh(alpha)*dbeta
    ddbddbeta = -2./tanh(alpha)*dalpha
    ddgddbeta = sin(2.*beta)*dgamma

    #------------

    daddgamma = 2./sinh(alpha)**2.*(dalpha*dgamma)
    dbddgamma = 2./sin(beta)**2.*(dbeta*dgamma)
    dgddgamma = 0.

    ddaddgamma = -2./tanh(alpha)*dgamma
    ddbddgamma = -2./tan(beta)*dgamma
    ddgddgamma = -2.*(dalpha/tanh(alpha) + dbeta/tan(beta))


    return np.array([
        [0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 1.],
        [daddalpha,dbddalpha,dgddalpha, ddaddalpha,ddbddalpha,ddgddalpha],
        [daddbeta ,dbddbeta ,dgddbeta , ddaddbeta ,ddbddbeta ,ddgddbeta ],
        [daddgamma,dbddgamma,dgddgamma, ddaddgamma,ddbddgamma,ddgddgamma]
    ])

def hypDistance(pos1, pos2):

    u = coords(pos1)
    v = coords(pos2)

    return hyperboloidDistance(u,v)

##### Distance Functions for coupling potentials and length constraints
def D12(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return np.cosh(a1)*np.cosh(a2) - np.sinh(a1)*np.cos(b1)*np.sinh(a2)*np.cos(b2) - np.sinh(a1)*np.sin(b1)*np.sinh(a2)*np.sin(b2)*np.cos(g1 - g2)

# First Derivatives of distance function
def da1D12(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return np.sinh(a1)*np.cosh(a2) - np.cosh(a1)*np.cos(b1)*np.sinh(a2)*np.cos(b2) - np.cosh(a1)*np.sin(b1)*np.sinh(a2)*np.sin(b2)*np.cos(g1 - g2)

def db1D12(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return np.sinh(a1)*np.sin(b1)*np.sinh(a2)*np.cos(b2) - np.sinh(a1)*np.cos(b1)*np.sinh(a2)*np.sin(b2)*np.cos(g1 - g2) 

def dg1D12(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return np.sinh(a1)*np.sin(b1)*np.sinh(a2)*np.sin(b2)*np.sin(g1 - g2)
# For the remaining three functions use:
# da2D12 = da1D12(state2, state1)
# db2D12 = db1D12(state2, state1)
# dg2D12 = dg1D12(state2, state1)

# Second Derivatives of distance function (Only needed for implicit methods for jacobian construction)
def da1D12a1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

def db1D12a1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return cosh(a1)*sin(b1)*sinh(a2)*cos(b2) - cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

def dg1D12a1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

def da2D12a1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return sinh(a1)*sinh(a2) - cosh(a1)*cos(b1)*cosh(a2)*cos(b2) - cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

def db2D12a1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return cosh(a1)*cos(b1)*sinh(a2)*sin(b2) - cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

def dg2D12a1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return -cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

def db1D12b1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return sinh(a1)*cos(b1)*sinh(a2)*cos(b2) + sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

def dg1D12b1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

def db2D12b1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

def dg2D12b1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return -sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

def dg1D12g1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

def dg2D12g1(state1, state2):
    a1,b1,g1 = state1.pos.copy()
    a2,b2,g2 = state2.pos.copy()
    return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)
# For the remaining nine functions of the upper triangular matrix use:
# da2D12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)

# da2D12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
# db2D12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)

# da2D12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
# db2D12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
# dg2D12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

# db2D12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
# dg2D12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

# dg2D12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)

# Expression to generate coupling potential terms in jacobian (i.e. second derivative of coupling potential)
def da2da1V12(m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1):
    # negative sign here so that the term can be added to geo terms later
    return -k/(m*f*sqrt( d12**2. - 1. ))*( (da1d12*da2d12)/sqrt( d12**2. - 1.) + ( arccosh(d12) - l )*( da2d12da1 - da1d12*(da2f/f + d12*da2d12/(d12**2. - 1.)) ) )


hypFuncDict = {
    "vertex_jac" : hyp_gen_jacobian,

    "d12"    : D12,

    "da1d12" : da1D12,
    "db1d12" : db1D12,
    "dg1d12" : dg1D12,

    "da1da1d12" : da1D12a1,
    "db1da1d12" : db1D12a1,
    "dg1da1d12" : dg1D12a1,
    "da2da1d12" : da2D12a1,
    "db2da1d12" : db2D12a1,
    "dg2da1d12" : dg2D12a1,

    "db1db1d12" : db1D12b1,
    "dg1db1d12" : dg1D12b1,
    "db2db1d12" : db2D12b1,
    "dg2db1d12" : dg2D12b1,

    "dg1dg1d12" : dg1D12g1,
    "dg2dg1d12" : dg2D12g1,

    "coupling_derivative" : da2da1V12
    }


hypSpace = Geometry(
    hypMetricTensor,
    hypChristoffel,
    hypDistance,
    hypFuncDict
    )

# // -------------------------------------------------------------
# // Model Information
# // -------------------------------------------------------------

def toPoincareBall(coord):

    p = coords(coord)

    x, y, z, w = p

    return np.array([x,y,z]) / (1. + w)

def pbScaling(pos):
    modelScale = 6.
    len2 = pos[0]**2. + pos[1]**2. + pos[2]**2.
    scale = modelScale**2. - len2
    return 4. * scale / (modelScale**2.)

hypModel = Model(
    toPoincareBall,
    pbScaling
    )

# // -------------------------------------------------------------
# // Obstacle/Bounding Ball Information
# // -------------------------------------------------------------

# // Sphere Bounding Box

# Default bounding box size (radius)
obstacleSize =2.

def distToSphere(pos):
    # //center point of H3 in coordinates:
    center = np.zeros(3)
    # //distance from center to position
    dist = hypDistance(pos.copy(),center)
    # //how far is this from the boundary sphere of radius 5?
    return obstacleSize - dist


sphereObstacle = Obstacle(
    distToSphere,
    obstacleSize
)

# //package stuff up for export
hyperbolic = AmbientSpace( hypSpace, hypModel, sphereObstacle)