# //the annoying thing I need to learn how to fix:
# //am importing a GLOBAL ambient space into the files defining my classes
# // :-(
from src.AmbientSpace import AmbientSpace

import numpy as np


# //right now the important parameters for the configuration space
# //are the list of masses and the radii of the balls involved
# //a state of the configuration space is a DataList of states of individual balls
# //a dState of the configuration space is also a DataList of dStates


class ConfigurationSpace:
    """
    A class used to store information about configuration space of system

    ...

    Attributes
    ----------
    masses : array
        array of masses for each vertex in the system

    radii : array
        array of radii of pseudosphere of each vertex in the system

    ambientspace : object
        AmbientSpace object describing ambient space containing system

    Methods
    -------
    dot(dataList1, dataList2)
        Calculate the dot product between dataList1 and dataList2 (statewise) with respect to the kinetic energy metric
        Returns scalar dot product value

    norm(dataList)
        Calculate the norm of each element in dataList with respect to kinetic energy metric
        Returns scalar norm value

    normalize(dataList)
        Calculate normalized dataList (normalize each state)
        Returns DataList object of normalized states from dataList

    obstacleCollisions(dataList)
        Generate list of vertices colliding with obstacles in the ambient space
        Returns array of indices cooresponding with which vertices have collided with obstacles

    obstacleGradient(dataList, indices)
        Compute the gradient of distance function to the boundary (vertex-obstacle collision) for only collided vertices in dataList listed by indices
        Returns DataList object clone of dataList with velocity of collided vertices updated via kinetic energy gradient

    ballCollisions(dataList)
        Generate list of vertices in dataList colliding with each other
        Returns array of indices cooresponding with which vertices have collided with each other

    ballGradient(dataList, indices)
        Compute the gradient of distance function to the boundary (vertex-vertex collision) for only collided vertices in dataList listed by indices
        Returns DataList object clone of dataList with velocity of collided vertices updated via kinetic energy gradient

    boundaryNormal(dataList, ball_collisionIndices, obstacle_collisionIndices)
        Calculate normal to boundary corresponding to all collision events (vertex-obstacle and vertex-vertex collision)
        Returns DataList object corresponding to all collision gradients (normalized)

    momentum(dataList)
        Calculate total momentum of all vertices in dataList (only for debugging in E3)
        Returns scalar of total momentum

    reflectIn(dataList, normal)
        Generate a clone of dataList updated corresponding to reflecting with boundary with normal vector given by normal
        Returns DataList object of reflecting dataList via normal to boundary
    """


    def __init__(self, masses, radii, ambientspace):
        self.N = len(masses)
        self.masses = masses
        self.radii = radii
        self.ambientspace = ambientspace


    # //get the dot product with respect to the kinetic energy metric
    def dot(self, dataList1, dataList2 ):
        dot = 0
        for i in range(self.N):
            # //add up the dot product for each individual ball
            dot += 1/2 * self.masses[i] * self.ambientspace.dot(dataList1.data[i],dataList2.data[i])
        return dot

    # //norm of the kinetic energy metric
    def norm(self, dataList ):
        return np.sqrt(self.dot(dataList, dataList))

    def normalize(self, dataList ):
        len = self.norm(dataList)
        res = dataList.clone().multiplyScalar(1/len)
        return res


    # //detect which balls collide with the ambient space's obstacle.
    # //output the indices of these balls as a list:
    def obstacleCollisions(self, dataList ):
        indices = []
        for i in range(self.N):
            posi = dataList.data[i].pos
            # print(posi)
            disti = self.ambientspace.distToObstacle( posi )
            # print("distance from wall")
            # print(disti)
            radi = self.radii[i]
            if( disti < radi ):
                # //the balls is intersecting the boundary:
                # //but, see if it is heading outward or inward
                newPos = dataList.data[i].clone().flow(0.001).pos
                newDist = self.ambientspace.distToObstacle(newPos)
                # //if this new distance is less, it's an intersection with inadmissable tangent
                if(newDist<disti):
                    indices.append(i)
        if(len(indices)==0):
            return indices
        print("obstacle collision")
        print(indices)
        return indices

    # //compute the gradient of the distance function to the boundary
    # //only at the specified balls: the rest are zero
    # //(instead; could take the gradient of the overall distance function which
    # //is a minimum over all distances to boundary, but this would have TONS
    # //of unnecessary computation)
    def obstacleGradient(self, dataList, indices):

        # //make a new state with same positions but zeroed out velocity:
        grad = dataList.clone()
        for i in range(self.N):
            grad.data[i].vel=np.array([0,0,0])

        # print(indices)
        if(len(indices) != 0):
            # //replace the velocity with the gradient in the correct index slots
            for index in range(len(indices)):
                i = indices[index]
                posi = dataList.data[i].pos

                # //with respect to the metric g on the ambient space X
                geomGradi = self.ambientspace.gradient(self.ambientspace.obstacle.distance, posi)

                # //the kinetic energy metric is 1/2*m*g, so the inverse metric tensor
                # //is scaled by 2/m:
                gradi = geomGradi.clone().multiplyScalar(2/self.masses[i])

                # //replace this in the gradient list:
                grad.data[i] = gradi

        return grad



    def ballCollisions(self, dataList ):
        indices = []

        for i in range(self.N):
            for j in range(self.N):
                if i != j:
                    distij = self.ambientspace.distance(dataList.data[i].pos, dataList.data[j].pos)
                    # print("ball to ball distance")
                    # print(distij)
                    radij = self.radii[i]+self.radii[j]

                    if(distij<radij):
                        # //the balls are intersecting: but are they approaching or moving apart?
                        newPosi = dataList.data[i].clone().flow(0.001).pos
                        newPosj = dataList.data[j].clone().flow(0.001).pos
                        newDist = self.ambientspace.distance(newPosi,newPosj)
                        # //if this new distance is less, it's an intersection with inadmissable tangent
                        if(newDist<distij):
                            indices.append([i,j])

        if( len(indices) == 0 ):
            return indices

        print("ball collision")
        print(indices)
        return indices

    # //compute the gradient of the distance function
    # //only compute for specified pairs [i,j] of balls, then add: rest are zero.
    def ballGradient(self, dataList, indices):

        # //make a new state with same positions but zeroed out velocity:
        grad = dataList.clone()
        for n in range(self.N):
            grad.data[n].vel=np.array([0,0,0])

        # //replace the velocity with the gradient in the correct index slots
        if(len(indices) != 0):
            for index in range(len(indices)):

                ij = indices[index]
                i = ij[0]
                j = ij[1]

                posi = dataList.data[i].pos
                posj = dataList.data[j].pos

                # //distance function to the ball "i":
                def disti(pos):
                    return self.ambientspace.distance(posi, pos)

                # //the gradient of this function, evaluated at position j
                gradjdisti = self.ambientspace.gradient(disti, posj)

                # //the kinetic energy metric for the jth particle is the Riemannian metric g,
                # //scaled by 1/2*m: thus the gradient is the g-gradient scaled by 2/m
                # //replace the gradient at j with this:
                grad.data[j] = gradjdisti.clone().multiplyScalar(2/self.masses[j])

                # //distance function to the ball "j":
                def distj(pos):
                    return self.ambientspace.distance(posj, pos)

                # //the gradient of this function, evaluated at position j
                gradidistj = self.ambientspace.gradient(distj, posi)

                # //the kinetic energy metric for the jth particle is the Riemannian metric g,
                # //scaled by 1/2*m: thus the gradient is the g-gradient scaled by 2/m
                # //replace the gradient at i with this:
                grad.data[i] = gradidistj.clone().multiplyScalar(2/self.masses[i])

        return grad


    def boundaryNormal(self, dataList, ball_collisionIndices, obstacle_collisionIndices):
        # print(obstacle_collisionIndices)
        grad1 = self.obstacleGradient(dataList, obstacle_collisionIndices)
        grad2 = self.ballGradient(dataList, ball_collisionIndices)

        # //add them together
        grad = grad1.clone().add(grad2)

        return self.normalize(grad.clone())

    # //this is MEANINGLESS outside of Euclidean space:
    # //this is only used for debugging: to confirm momentum is conserved
    # //with ball-on-ball collisions (but not collisions with the boundary, obv)
    def momentum(self, dataList):
        p = np.array([0,0,0])
        for i in range(self.N):
            p.add(dataList.data[i].vel.clone().multiplyScalar(self.masses[i]))

        return p


    # //reflect a state in a normal vector
    # //dataList is the current tangent vector to configuration space (dataList of all the balls)
    # //normal is the normal vector to the boundary of configuration space
    def reflectIn(self, dataList, normal):

        dot = self.dot(dataList,normal)
        norm2 = self.dot(normal,normal)

        coef = 2.*dot/norm2

        result =  dataList.clone().sub(normal.clone().multiplyScalar(coef))

        return result
