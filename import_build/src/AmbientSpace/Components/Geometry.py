import numpy as np

from src.Computation.State import State
from src.Computation.dState import dState


# Class for storing ambient space geometry
class Geometry:

    def __init__(self, metricTensor, christoffel, distance):
        self.metricTensor = metricTensor
        self.christoffel = christoffel
        self.distance = distance


    # Take covariant derivative of a state
    # return dState storing the initial velocity, and the covariant coordinates (x'', y'', z'')= of the covariant derivative
    def covariantAcceleration(self, state):
        vel = state.vel
        acc = self.christoffel(state)
        return dState(vel,acc)

    # Calculate the dot product of two vectors, using the  metric tensor
    def dot(self, state1, state2):
        mat = self.metricTensor(state1.pos)

        v1 = state1.vel.copy()
        v2 = state2.vel.copy()

        # Compute the dot product
        return v1 @ (v2 @ mat)


    # Get a basis for the tangent space at a point ( in coordinates )
    # @@@@ Consider changing how we want to do this, depending on which implementation of the gradient we want below:
    # @@@@ right now, this is the coordinate basis, and we use the metric tensor
    # @@@@ in the gradient: could instead do Gram-Schmidt here then calculate
    # @@@@ gradient as differential directly in that basis.
    def tangentBasis(pos):
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
        #//using the hyperbolic metric:
        metric = self.metricTensor(pos)
        if(abs(metric.determinant())<0.00001):
            print('Warning! Metric Tensor Near Singular')
            print(pos)
            print(metric)

        invMetric = metric.clone().invert()
        differential.vel = differential.vel.clone().applyMatrix3(invMetric)

        return differential