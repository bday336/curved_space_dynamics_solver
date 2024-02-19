####################
# Helper Functions #
####################

import numpy as np
from numpy import zeros,array,arange,sqrt,sin,cos,tan,sinh,cosh,tanh,pi,arcsinh,arccosh,arctanh,arccos,arctan2,matmul,exp,identity,append

######################
# H3 Space Functions #
######################

# SO(3,1) -> For H3 isometries (x,y,z,w)

def boostxh3(u): 
	return array([
	   [cosh(u), 0., 0., sinh(u)],
	   [0., 1., 0., 0.],
	   [0., 0., 1., 0.],
	   [sinh(u), 0., 0., cosh(u)]
		])

def boostyh3(u): 
	return array([
	   [1, 0., 0., 0],
	   [0., cosh(u), 0., sinh(u)],
	   [0., 0., 1., 0.],
	   [0., sinh(u), 0., cosh(u)]
		])

def rotxh3(u): 
	return array([
	   [1., 0., 0., 0.],
	   [0., cos(u), -sin(u), 0.],
	   [0., sin(u), cos(u), 0.],
       [0., 0., 0., 1.]
		])

def rotyh3(u): 
	return array([
	   [cos(u), 0., sin(u), 0.],
	   [0., 1., 0., 0.],
	   [-sin(u), 0., cos(u), 0.],
       [0., 0., 0., 1.]
		])      

def rotzh3(u): 
	return array([
	   [cos(u), -sin(u), 0., 0.],
	   [sin(u), cos(u), 0., 0.],
       [0., 0., 1., 0.],
       [0., 0., 0., 1.]
		])

def rot2transh3pos(rot_vec):
    # Coordinate Killing axis along x-axis
    al,be,ga = rot_vec
    a = arcsinh(sinh(al)*sin(be)*sin(ga))
    b = arcsinh(sinh(al)*sin(be)*cos(ga)/cosh(a))
    g = arctanh(tanh(al)*cos(be))
    return array([a,b,g])

def rot2transh3vel(rot_vec_p,rot_vec_v):
    # Coordinate Killing axis along x-axis
    al,be,ga = rot_vec_p
    dal,dbe,dga = rot_vec_v
    a,b,g = rot2transh3pos(rot_vec_p)
    da = (dal*cosh(al)*sin(be)*sin(ga) + dbe*sinh(al)*cos(be)*sin(ga) + dga*sinh(al)*sin(be)*cos(ga))/cosh(a)
    db = (dal*cosh(al)*sin(be)*cos(ga) + dbe*sinh(al)*cos(be)*cos(ga) - dga*sinh(al)*sin(be)*sin(ga) - da*sinh(a)*sinh(b))/(cosh(a)*cosh(b))
    dg = (dal*(cosh(al)*cos(be)*cosh(g) - sinh(al)*sinh(g)) - dbe*sinh(al)*sin(be)*cosh(g))/(cosh(a)*cosh(b))
    return array([da,db,dg])


def rot2hyp(vec):
    a,b,g = vec
    return np.array([
        sinh(a)*sin(b)*cos(g), 
        sinh(a)*sin(b)*sin(g), 
        sinh(a)*cos(b), 
        cosh(a)])

def hyp2rot(vec): 
    x,y,z,w = vec
    return np.array([
        arccosh(w), 
        arccos(z/sinh(arccosh(w))), 
        arctan2(y, x)])

def hyp2poin3d(vec):
    x,y,z,w = vec
    return np.array([
        x/(w + 1), 
        y/(w + 1),
        z/(w + 1)])

def h3dist(v1,v2):
    x1,y1,z1,w1 = v1
    x2,y2,z2,w2 = v2
    return arccosh(-x1*x2 - y1*y2 - z1*z2 + w1*w2)

def h3distb2(v1,v2):
    x1,y1,z1,w1 = v1
    x2,y2,z2,w2 = v2
    return -x1*x2 - y1*y2 - z1*z2 + w1*w2

def killingvech3(pos, v, dir):
    x,y,z,w = rot2hyp(pos)
    if(dir == "x"):
        killingv = np.array([v*w, 0., 0., v*x])
    elif(dir == "y"):
        killingv = np.array([0., v*w, 0., v*y])
    elif(dir == "z"):
        killingv = np.array([0., 0., v*w, v*z])
    elif(dir == "vz"):
        killingv = np.array([v*y, -v*x, 0., 0.])

    ad = killingv[3]/sinh(pos[0])
    bd = (ad*cosh(pos[0])*cos(pos[1]) - killingv[2])/(sinh(pos[0])*sin(pos[1]))
    gd = (killingv[1]*cos(pos[2]) - killingv[0]*sin(pos[2]))/(sinh(pos[0])*sin(pos[1]))
    return np.array([ad, bd, gd])

# Use Rotational Parameterization
def hypercirch3(center,rad):
    u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
    x = sinh(rad)*sin(u)*cos(v)
    y = sinh(rad)*sin(u)*sin(v)
    z = sinh(rad)*cos(u)
    w = np.full(u.shape,cosh(rad))

    for b in range(u.shape[1]):
        for a in range(u.shape[0]):
            # This method currently does not preserve the orientation of the sphere (need to update if we wanted to have some applied texture)
            testarr=rotxh3(arctan2(center[2],center[1])) @ rotzh3(arctan2((rotxh3(-arctan2(center[2],center[1])) @ center)[1],(rotxh3(-arctan2(center[2],center[1])) @ center)[0])) @ boostxh3(arccosh(center[3])) @ array([x[b,a],y[b,a],z[b,a],w[b,a]])
            x[b,a],y[b,a],z[b,a],w[b,a]=testarr[0],testarr[1],testarr[2],testarr[3]

    return x/(w+1.),y/(w+1.),z/(w+1.)


######################
# S3 Space Functions #
######################

def rot2r4(vec): 
    a,b,g = vec
    return np.array([
        sin(a)*sin(b)*cos(g),
        sin(a)*sin(b)*sin(g), 
        sin(a)*cos(b),
        cos(a)])

def r42rot(vec): 
    x,y,z,w = vec
    return np.array([
        arccos(w), 
        arccos(z/sin(arccos(w))),
        arctan2(y, x)])

def s2rstproj(vec):
    x,y,z,w = vec
    return np.array([
        x/(1. - w),
        y/(1. - w), 
        z/(1. - w)])

def r22s2stproj(vec): 
    xp,yp = vec
    return np.array([
        (2.*xp)/(1. + xp**2. + yp**2.),
        (2.*yp)/(1. + xp**2. + yp**2.),
        (xp**2. + yp**2. - 1.)/(1. + xp**2. + yp**2.)])

def r4dist(v1, v2):
    x1,y1,z1,w1 = v1
    x2,y2,z2,w2 = v2
    return arccos(x1*x2 + y1*y2 + z1*z2 + w1*w2)

def r4distb2(v1, v2):
    x1,y1,z1,w1 = v1
    x2,y2,z2,w2 = v2
    return x1*x2 + y1*y2 + z1*z2 + w1*w2

def killingvecs3(pos, v, dir):
    x,y,z,w = rot2r4(pos)
    if(dir == "vx"):
        killingv = np.array([0., v*z, -v*y, 0.])
    elif(dir == "vy"):
        killingv = np.array([-v*z, 0., v*x, 0.])
    elif(dir == "vz"):
        killingv = np.array([v*y, -v*x, 0., 0.])
    elif(dir == "bx"):
        killingv = np.array([v*w, 0., 0., -v*x])
    elif(dir == "by"):
        killingv = np.array([0., -v*w, 0., v*y])
    elif(dir == "bz"):
        killingv = np.array([0., 0., v*w, -v*z])
    ad = -killingv[3]/sin(pos[0])
    bd = (ad*cos(pos[0])*cos(pos[1]) - killingv[2])/(sin(pos[0])*sin(pos[1]))
    gd = (killingv[1]*cos(pos[2]) - killingv[0]*sin(pos[2]))/(sin(pos[0])*sin(pos[1]))
    return np.array([ad, bd, gd])


##########################
# Import Build Functions #
##########################

def genballe3(center,rad):
    u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
    x = rad*np.sin(u)*np.cos(v)
    y = rad*np.sin(u)*np.sin(v)
    z = rad*np.cos(u)

    for b in range(u.shape[1]):
        for a in range(u.shape[0]):
            # This method currently does not preserve the orientation of the sphere (need to update if we wanted to have some applied texture)
            testarr=np.array([x[b,a],y[b,a],z[b,a]]) + center
            x[b,a],y[b,a],z[b,a]=testarr[0],testarr[1],testarr[2]

    return x,y,z

def genballh3(center,rad):
    u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
    x = np.sinh(rad)*np.sin(u)*np.cos(v)
    y = np.sinh(rad)*np.sin(u)*np.sin(v)
    z = np.sinh(rad)*np.cos(u)
    w = np.full(u.shape,np.cosh(rad))

    for b in range(u.shape[1]):
        for a in range(u.shape[0]):
            # This method currently does not preserve the orientation of the sphere (need to update if we wanted to have some applied texture)
            testarr=rotxh3(np.arctan2(center[2],center[1])) @ rotzh3(np.arctan2((rotxh3(-np.arctan2(center[2],center[1])) @ center)[1],(rotxh3(-np.arctan2(center[2],center[1])) @ center)[0])) @ boostxh3(np.arccosh(center[3])) @ np.array([x[b,a],y[b,a],z[b,a],w[b,a]])
            x[b,a],y[b,a],z[b,a],w[b,a]=testarr[0],testarr[1],testarr[2],testarr[3]

    return x/(w+1.),y/(w+1.),z/(w+1.)

