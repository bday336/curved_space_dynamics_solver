import numpy as np
import random

def randomVec3Sphere( Radius ):
    theta = random.uniform(0,6.29)
    z = random.uniform(-1,1)

    pt = np.array( [np.cos(theta), np.sin(theta), z] ) * Radius

    return pt


def randomVec3Ball( Radius ):
    pt = randomVec3Sphere( Radius )

    r = random.uniform(0,1)
    r = r**1.3333

    return (pt * r)

