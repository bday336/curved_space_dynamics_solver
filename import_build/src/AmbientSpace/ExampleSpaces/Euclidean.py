import numpy as np

from src.AmbientSpace.Components.Geometry import Geometry
from src.AmbientSpace.Components.Model import Model
from src.AmbientSpace.Components.Obstacle import Obstacle

from src.AmbientSpace.AmbientSpace import AmbientSpace

from src.Computation.State import State

# // -------------------------------------------------------------
# //some euclidean geometry stuff: (default to 3D)
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
# //model of Euclidean space : do nothing
# // -------------------------------------------------------------


def identityR3(coords):
    return coords

def unitScaling(pos):
    return 1.


eucModel = Model(identityR3,unitScaling)





# // -------------------------------------------------------------
# //obstacles to contain balls in Euclidean Space
# // -------------------------------------------------------------

# //a sphere
R = 6.

def distToSphere(pos):
    return R-np.sqrt(pos[0]**2. + pos[1]**2. + pos[2]**2.)

sphereGeom = None #new SphereBufferGeometry(R,64,32);

generateSphState = None #function(){
#     let pos = randomVec3Ball(0.8*R);
#     # let vel = randomVec3Ball(1);
#     return new State(pos,vel);
# }


sphereObstacle = Obstacle(
    distToSphere,
    sphereGeom,
    R,
    generateSphState
)







# //a box

def distToBox(pos):
    xWall = 6. - abs(pos.x)
    yWall = 4. - abs(pos.y)
    zWall = 4. - abs(pos.z)

    return min(xWall,min(yWall,zWall))


boxGeom = None #new BoxBufferGeometry(12,8,8);

boxSize = 6.

generateBoxState = None #function(){
#     let x = 6.*Math.random()-6;
#     let y = 4.*Math.random()-4;
#     let z = 4.*Math.random()-4;
#     let pos = new Vector3(x,y,z).multiplyScalar(0.8);
#     let vel = randomVec3Ball(1);
#     return new State(pos,vel);
# }



boxObstacle = Obstacle(
    distToBox,
    boxGeom,
    boxSize,
    generateBoxState,
)






# //package stuff up for export
euclidean = AmbientSpace( eucSpace, eucModel, sphereObstacle)

