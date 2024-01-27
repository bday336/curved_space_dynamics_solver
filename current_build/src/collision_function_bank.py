#######################
# Collision Functions #
#######################

# Collision functions for two rigid vertices

import numpy as np
from numpy import zeros,array,arange,sqrt,sin,cos,tan,sinh,cosh,tanh,pi,arcsinh,arccosh,arctanh,arccos,arctan2,matmul,exp,identity,append
from function_bank import rot2hyp, rot2r4, boostxh3, rotxh3, rotyh3, rotzh3

######################
# H3 Space Functions #
######################

def h3collision(state_vec,params):
    def D12(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)
    
    # First Derivatives
    def da1D12(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*cosh(a2) - cosh(a1)*cos(b1)*sinh(a2)*cos(b2) - cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def db1D12(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*sin(b1)*sinh(a2)*cos(b2) - sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2) 

    def dg1D12(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)
    # For the remaining three functions use:
    # da2D12 = da1D12(a2, b2, g2, a1, b1, g1)
    # db2D12 = db1D12(a2, b2, g2, a1, b1, g1)
    # dg2D12 = dg1D12(a2, b2, g2, a1, b1, g1)

    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    m1,m2,r1,r2 = params


    # Distance Function
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # First derivatives of distance function
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    print(d12,da1d12,db1d12,dg1d12,da2d12,db2d12,dg2d12)

    da1 = (2*da1d12)/(m1*sqrt(d12**2. - 1.))
    db1 = (2*db1d12)/(m1*sinh(a1)**2.*sqrt(d12**2. - 1.))
    dg1 = (2*dg1d12)/(m1*sinh(a1)**2.*sin(b1)**2.*sqrt(d12**2. - 1.))
    da2 = (2*da2d12)/(m2*sqrt(d12**2. - 1.))
    db2 = (2*db2d12)/(m2*sinh(a2)**2.*sqrt(d12**2. - 1.))
    dg2 = (2*dg2d12)/(m2*sinh(a2)**2.*sin(b2)**2.*sqrt(d12**2. - 1.))

    vels = np.array([ad1, bd1, gd1, ad2, bd2, gd2])

    grad = np.array([da1, db1, dg1, da2, db2, dg2])

    gradnorm = np.sqrt(m1/2.*da1**2. + m1*sinh(a1)**2./2.*db1**2. + m1*sinh(a1)**2.*sin(b1)**2./2.*dg1**2. + 
                       m2/2.*da2**2. + m2*sinh(a2)**2./2.*db2**2. + m2*sinh(a2)**2.*sin(b2)**2./2.*dg2**2.)

    reflection = vels - (2./(gradnorm**2.))*(m1/2.*ad1*da1 + m1*sinh(a1)**2./2.*bd1*db1 + m1*sinh(a1)**2.*sin(b1)**2./2.*gd1*dg1 + m2/2.*ad2*da2 + m2*sinh(a2)**2./2.*bd2*db2 + m2*sinh(a2)**2.*sin(b2)**2./2.*gd2*dg2) * grad

    return np.concatenate((np.array([a1,b1,g1,a2,b2,g2]),reflection))

def rot2hypv(pos, vel): 
    return np.array([
        vel[0]*cosh(pos[0])*sin(pos[1])*cos(pos[2]) + vel[1]*sinh(pos[0])*cos(pos[1])*cos(pos[2]) - vel[2]*sinh(pos[0])*sin(pos[1])*sin(pos[2]),
        vel[0]*cosh(pos[0])*sin(pos[1])*sin(pos[2]) + vel[1]*sinh(pos[0])*cos(pos[1])*sin(pos[2]) + vel[2]*sinh(pos[0])*sin(pos[1])*cos(pos[2]),
        vel[0]*cosh(pos[0])*cos(pos[1]) - vel[1]*sinh(pos[0])*sin(pos[1]),
        vel[0]*sinh(pos[0])])

def hyp2rotv(pos, vel): 
    return np.array([
        vel[3]/sinh(pos[0]),
        (vel[3]/tanh(pos[0])*cos(pos[1]) - vel[2])/(sinh(pos[0])*sin(pos[1])),
        (vel[1]*cos(pos[2]) - vel[0]*sin(pos[2]))/(sinh(pos[0])*sin(pos[1]))])

def h3kvecproj(pos, v, dir): 
    hyppos = rot2hyp(pos)
    hypvel = rot2hypv(pos, v)
    if (dir == "bx"):
        kvec = np.array([hyppos[3], 0., 0., hyppos[0]])
    if (dir == "by"):
        kvec = np.array([0., hyppos[3], 0., hyppos[1]])
    if (dir == "bz"):
        kvec = np.array([0., 0., hyppos[3], hyppos[2]])
    if (dir == "wx"):
        kvec = np.array([0., -hyppos[2], hyppos[1], 0.])
    if (dir == "wy"):
        kvec = np.array([hyppos[2], 0., -hyppos[0], 0.])
    if (dir == "wz"):
        kvec = np.array([-hyppos[1], hyppos[0], 0., 0.])
    # print(-kvec[0]**2. - kvec[1]**2. - kvec[2]**2. + kvec[3]**2.)
    kvecnorm = np.sqrt(abs(kvec[0]**2. + kvec[1]**2. + kvec[2]**2. - kvec[3]**2.))
    kproj = (hypvel[0]*kvec[0] + hypvel[1]*kvec[1] + hypvel[2]*kvec[2] - hypvel[3]*kvec[3])/kvecnorm
    kprojvec = kproj*kvec/kvecnorm
    return kproj,kprojvec

def transform2origh3col(pos1, pos2, vel1, vel2, dist): 
    pos1hyp = rot2hyp(pos1)
    pos2hyp = rot2hyp(pos2)
    vel1hyp = vel1
    vel2hyp = vel2

    # Assume in xy plane
    transform1 = rotzh3(arctan2(pos1hyp[1], pos1hyp[0])) @ boostxh3(-arccosh(pos1hyp[3])) @ rotzh3(-arctan2(pos1hyp[1], pos1hyp[0]))
    trans12op1 = transform1 @ pos1hyp
    trans12op2 = transform1 @ pos2hyp
    trans12ov1 = transform1 @ vel1hyp
    trans12ov2 = transform1 @ vel2hyp
    print('transform 1')
    print(trans12op1)
    print(trans12op2)
    print(trans12ov1)
    print(trans12ov2)
    
    # transform2 = rotzh3(arctan2(pos1hyp[1], pos1hyp[0])) @ rotxh3(np.pi/2.) @ rotzh3(arctan2(trans12xyp1[1], trans12xyp1[0])) @ boostxh3(-arccosh(trans12xyp1[3])) @ rotzh3(-arctan2(trans12xyp1[1], trans12xyp1[0]))
    # trans12op1 = transform2 @ trans12xyp1
    # trans12op2 = transform2 @ trans12xyp2
    # trans12ov1 = transform2 @ trans12xyv1
    # trans12ov2 = transform2 @ trans12xyv2
    # print('transform 2')
    # print(trans12op1)
    # print(trans12op2)
    # print(trans12ov1)
    # print(trans12ov2)
    
    transform3 = rotzh3(-arctan2(trans12op2[1], trans12op2[0]))
    trans22xp1 = transform3 @ trans12op1
    trans22xp2 = transform3 @ trans12op2
    trans22xv1 = transform3 @ trans12ov1
    trans22xv2 = transform3 @ trans12ov2
    print('transform 3')
    print(trans22xp1)
    print(trans22xp2)
    print(trans22xv1)
    print(trans22xv2)

    # transform4 = rotzh3(-arctan2(trans22xyp2[1], trans22xyp2[0]))
    # trans22xp1 = transform4 @ trans22xyp1
    # trans22xp2 = transform4 @ trans22xyp2
    # trans22xv1 = transform4 @ trans22xyv1
    # trans22xv2 = transform4 @ trans22xyv2
    # print('transform 4')
    # print(trans22xp1)
    # print(trans22xp2)
    # print(trans22xv1)
    # print(trans22xv2)

    transform5 = boostxh3(-.5*dist)
    transm2op1 = transform5 @ trans22xp1
    transm2op2 = transform5 @ trans22xp2
    transm2ov1 = transform5 @ trans22xv1
    transm2ov2 = transform5 @ trans22xv2
    print('transform 5')
    print(transm2op1)
    print(transm2op2)
    print(transm2ov1)
    print(transm2ov2)

    transv2ov1 = boostxh3(arccosh(transm2op1[3])) @ transm2ov1
    transv2ov2 = boostxh3(-arccosh(transm2op2[3])) @ transm2ov2
    print('at origin')
    print(transv2ov1)
    print(transv2ov2)

    totaldiff = transv2ov1 + transv2ov2

    # totaldiffcol = np.linalg.inv(transform)

    return transv2ov1,transv2ov2,totaldiff

######################
# S3 Space Functions #
######################

def s3collision(state_vec,params):
    def D12(a1, b1, g1, a2, b2, g2):
        return cos(a1)*cos(a2) + sin(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))
    
    # First Derivatives
    def da1D12(a1, b1, g1, a2, b2, g2):
        return -cos(a2)*sin(a1) + cos(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

    def db1D12(a1, b1, g1, a2, b2, g2):
        return sin(a1)*sin(a2)*(-cos(b2)*sin(b1) + cos(b1)*cos(g1 - g2)*sin(b2))

    def dg1D12(a1, b1, g1, a2, b2, g2):
        return -sin(a1)*sin(a2)*sin(b1)*sin(b2)*sin(g1 - g2)
    # For the remaining three functions use:
    # da2D12 = da1D12(a2, b2, g2, a1, b1, g1)
    # db2D12 = db1D12(a2, b2, g2, a1, b1, g1)
    # dg2D12 = dg1D12(a2, b2, g2, a1, b1, g1)
    
    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    m1,m2,r1,r2 = params

    # Distance Function
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # First derivatives of distance function
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)

    da1 = -(2.*da1d12)/(m1*sqrt(1. - d12**2.))
    db1 = -(2.*db1d12)/(m1*sin(a1)**2. * sqrt(1. - d12**2.))
    dg1 = -(2.*dg1d12)/(m1*sin(a1)**2. * sin(b1)**2. * sqrt(1. - d12**2.))
    da2 = -(2.*da2d12)/(m2*sqrt(1. - d12**2.))
    db2 = -(2.*db2d12)/(m2*sin(a2)**2. * sqrt(1. - d12**2.))
    dg2 = -(2.*dg2d12)/(m2*sin(a2)**2. * sin(b2)**2. * sqrt(1. - d12**2.))

    vels = np.array([ad1, bd1, gd1, ad2, bd2, gd2])

    grad = np.array([da1, db1, dg1, da2, db2, dg2])

    gradnorm = np.sqrt(m1/2.*da1**2. + m1*sin(a1)**2./2.*db1**2. + m1*sin(a1)**2.*sin(b1)**2./2.*dg1**2. + 
                       m2/2.*da2**2. + m2*sin(a2)**2./2.*db2**2. + m2*sin(a2)**2.*sin(b2)**2./2.*dg2**2.)

    reflection = vels - (2./(gradnorm)**2.)*(m1/2.*ad1*da1 + m1*sin(a1)**2./2.*bd1*db1 + m1*sin(a1)**2.*sin(b1)**2./2.*gd1*dg1 + m2/2.*ad2*da2 + m2*sin(a2)**2./2.*bd2*db2 + m2*sin(a2)**2.*sin(b2)**2./2.*gd2*dg2) * grad

    return np.concatenate((np.array([a1,b1,g1,a2,b2,g2]),reflection))

def rot2r4v(pos, vel): 
    return np.array([
        vel[0]*cos(pos[0])*sin(pos[1])*cos(pos[2]) + vel[1]*sin(pos[0])*cos(pos[1])*cos(pos[2]) - vel[2]*sin(pos[0])*sin(pos[1])*sin(pos[2]),
        vel[0]*cos(pos[0])*sin(pos[1])*sin(pos[2]) + vel[1]*sin(pos[0])*cos(pos[1])*sin(pos[2]) + vel[2]*sin(pos[0])*sin(pos[1])*cos(pos[2]),
        vel[0]*cos(pos[0])*cos(pos[1]) - vel[1]*sin(pos[0])*sin(pos[1]),
        -vel[0]*sin(pos[0])])

def s3kvecproj(pos, v, dir): 
    r4pos = rot2r4(pos)
    r4vel = rot2r4v(pos, v)
    if(dir == "bx"):
        kvec = np.array([r4pos[3], 0., 0., -r4pos[0]])
    if(dir == "by"):
        kvec = np.array([0., r4pos[3], 0., -r4pos[1]])
    if(dir == "bz"):
        kvec = np.array([0., 0., r4pos[3], -r4pos[2]])
    if(dir == "vx"):
        kvec = np.array([0., -r4pos[2], r4pos[1], 0.])
    if(dir == "vy"):
        kvec = np.array([r4pos[2], 0., -r4pos[0], 0.])
    if(dir == "vz"):
        kvec = np.array([-r4pos[1], r4pos[0], 0., 0.])
    # print(-kvec[0]**2. - kvec[1]**2. - kvec[2]**2. + kvec[3]**2.)
    kvecnorm = np.sqrt(abs(kvec[0]**2. + kvec[1]**2. + kvec[2]**2. + kvec[3]**2.))
    kproj = (r4vel[0]*kvec[0] + r4vel[1]*kvec[1] + r4vel[2]*kvec[2] + r4vel[3]*kvec[3])/kvecnorm
    kprojvec = kproj*kvec/kvecnorm
    return kproj,kprojvec

######################
# E3 Space Functions #
######################

def e3collision(state_vec,params):
    def D12(a1, b1, g1, a2, b2, g2):
        return (a2-a1)**2. + (b2-b1)**2. + (g2-g1)**2.
    
    # First Derivatives
    def da1D12(a1, b1, g1, a2, b2, g2):
        return -2.*(a2-a1)

    def db1D12(a1, b1, g1, a2, b2, g2):
        return -2.*(b2-b1)

    def dg1D12(a1, b1, g1, a2, b2, g2):
        return -2.*(g2-g1)
    # For the remaining three functions use:
    # da2D12 = da1D12(a2, b2, g2, a1, b1, g1)
    # db2D12 = db1D12(a2, b2, g2, a1, b1, g1)
    # dg2D12 = dg1D12(a2, b2, g2, a1, b1, g1)

    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    m1,m2,r1,r2 = params


    # Distance Function
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # First derivatives of distance function
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)

    da1 = (sqrt(2)*da1d12)/(sqrt(m1)*2*sqrt(d12))
    db1 = (sqrt(2)*db1d12)/(sqrt(m1)*2*sqrt(d12))
    dg1 = (sqrt(2)*dg1d12)/(sqrt(m1)*2*sqrt(d12))
    da2 = (sqrt(2)*da2d12)/(sqrt(m2)*2*sqrt(d12))
    db2 = (sqrt(2)*db2d12)/(sqrt(m2)*2*sqrt(d12))
    dg2 = (sqrt(2)*dg2d12)/(sqrt(m2)*2*sqrt(d12))

    vels = np.array([ad1, bd1, gd1, ad2, bd2, gd2])

    grad = np.array([da1, db1, dg1, da2, db2, dg2])

    gradnorm = np.sqrt(m1/2.*da1**2. + m1/2.*db1**2. + m1/2.*dg1**2. + 
                       m2/2.*da2**2. + m2/2.*db2**2. + m2/2.*dg2**2.)

    reflection = vels - (2./(gradnorm**2.))*(m1/2.*ad1*da1 + m1/2.*bd1*db1 + m1/2.*gd1*dg1 + m2/2.*ad2*da2 + m2/2.*bd2*db2 + m2/2.*gd2*dg2) * grad

    return np.concatenate((np.array([a1,b1,g1,a2,b2,g2]),reflection))


def e3kvecproj(pos, v, dir): 
    if (dir == "x"):
        kvec = np.array([1, 0., 0.])
    if (dir == "y"):
        kvec = np.array([0., 1, 0.])
    if (dir == "z"):
        kvec = np.array([0., 0., 1])
    # print(-kvec[0]**2. - kvec[1]**2. - kvec[2]**2. + kvec[3]**2.)
    kvecnorm = np.sqrt(abs(kvec[0]**2. + kvec[1]**2. + kvec[2]**2.))
    kproj = (v[0]*kvec[0] + v[1]*kvec[1] + v[2]*kvec[2])/kvecnorm
    return kproj
