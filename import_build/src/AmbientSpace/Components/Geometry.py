import numpy as np

from src.Computation.State import State
from src.Computation.dState import dState


# A class used to store information about geometry of space
class Geometry:
    """
    A class used to store information about geometry of space

    ...

    Attributes
    ----------
    metricTensor : function
        function to generate 2D matrix representation of metric tensor at given position

    christoffel : function
        function to generate the terms corresponding to christoffel terms for acceleration of each vertex

    distance : function
        function to calculate the distance between two points in the ambient space

    Methods
    -------
    covariantAcceleration(state)
        Takes covariant derivative of state
        Returns dState storing initial velocity and covariant coordinates (second derivatives of coordinates) of the covariant derivatives

    dot(state1, state2)
        Calculates the dot product of velocity vectors for specificed states using the metric tensor
        Returns scalar value of dot product

    tangentBasis(pos)
        Calculates a basis for the tangent space at a point (in coordinates)
        Returns array of basis vectors

    gradient(fn, pos)
        Calculate the gradient of function fn at position pos
        Returns differential (gradient) of function at given point
    """

    def __init__(self, metricTensor, christoffel, distance):
        self.metricTensor = metricTensor
        self.christoffel  = christoffel
        self.distance     = distance


    def covariantAcceleration(self, state):
        vel = state.vel
        acc = self.christoffel(state)
        return dState(vel,acc)


    def dot(self, state1, state2):
        mat = self.metricTensor(state1.pos)
        v1 = state1.vel.copy()
        v2 = state2.vel.copy()
        # Apply this to the second vector
        gv2 = mat @ v2
        # Compute the dot product
        return v1.dot(gv2)


    # Get a basis for the tangent space at a point ( in coordinates )
    # @@@@ Consider changing how we want to do this, depending on which implementation of the gradient we want below:
    # @@@@ right now, this is the coordinate basis, and we use the metric tensor
    # @@@@ in the gradient: could instead do Gram-Schmidt here then calculate
    # @@@@ gradient as differential directly in that basis.
    def tangentBasis(self,pos):
        b1 = State(pos, np.array([1,0,0]))
        b2 = State(pos, np.array([0,1,0]))
        b3 = State(pos, np.array([0,0,1]))
        return np.array([b1,b2,b3])


    # //WARNING: IF THE COORDINATE SYSTEM IS SINGULAR: THIS COMPUTATION IS BAD AT THAT POINT!
    # //NEED GOOD COORDINATES.....
    # //get the gradient of a function fn at a position pos
    def gradient(self, fn, pos):

        basis = self.tangentBasis(pos)
        differential = State(pos, np.array([0,0,0]))

        #//add them all up:
        df0 = basis[0].differentiate(fn)
        b0 = basis[0].clone().multiplyScalar(df0)
        differential.add(b0)

        df1 = basis[1].differentiate(fn)
        b1 = basis[1].clone().multiplyScalar(df1)
        differential.add(b1)

        df2 = basis[2].differentiate(fn)
        b2 = basis[2].clone().multiplyScalar(df2)
        differential.add(b2)

        #//now the differential needs to be converted from a covector to a vector
        #//using the inverse metric:
        metric = self.metricTensor(pos)
        if(abs(np.linalg.det(metric))<0.00001):
            print('Warning! Metric Tensor Near Singular')
            print(pos)
            print(metric)

        invMetric = np.linalg.inv(metric.copy())
        differential.vel = invMetric @ differential.vel.copy()

        return differential