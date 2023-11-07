##########################
# System Setup Functions #
##########################

import numpy as np
from numpy import zeros,array,arange,sqrt,sin,cos,tan,sinh,cosh,tanh,pi,arccos,arcsinh,arccosh,arctanh,arctan2,matmul,exp,identity,append,linalg,add

#########################################################
#################### H3 Test Systems ####################
#########################################################

###################
# H3 Sim Geodesic #
###################

def dynfunc_h3simgeo(state_vec,params):
    
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,m = params


    riga1 = 1./2.*sinh(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) + 0.
    rigb1 = -2. * ad1*bd1/tanh(a1) + .5*sin(2.*b1)*gd1**2         + 0.
    rigg1 = -2. * ad1*gd1/tanh(a1) - 2.*bd1*gd1/tan(b1)           + 0.

    return np.array(
        [ad1,bd1,gd1,
        riga1,rigb1,rigg1])

# Setup system jacobian (Needed since method is implicit)
def dynjac_h3simgeo(state_vec, params):
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,m = params

    #---------- particle 1

    da1riga1 = cosh(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + 0.
    db1riga1 = cos(b1)*sin(b1)*sinh(2.*a1)*gd1**2.           + 0.
    dg1riga1 = 0.                                            + 0.

    dad1riga1 = 0.                          + 0.
    dbd1riga1 = sinh(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sinh(2.*a1)*gd1 + 0.

    #----------

    da1rigb1 = 2./sinh(a1)**2.*(ad1*bd1) + 0.
    db1rigb1 = cos(2.*b1)*gd1**2.        + 0.
    dg1rigb1 = 0.                        + 0.

    dad1rigb1 = -2./tanh(a1)*bd1 + 0.
    dbd1rigb1 = -2./tanh(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1    + 0.

    #------------

    da1rigg1 = 2./sinh(a1)**2.*(ad1*gd1) + 0.
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)  + 0.
    dg1rigg1 = 0.                        + 0.

    dad1rigg1 = -2./tanh(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                  + 0.
    dgd1rigg1 = -2.*(ad1/tanh(a1) + bd1/tan(b1)) + 0.


    return np.array([
        [0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, dad1riga1,dbd1riga1,dgd1riga1],
        [da1rigb1,db1rigb1,dg1rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1]
    ])



#######################
# H3 Exact Spring Bar #
#######################

def dynfunc_h3exactbar(state_vec,params):
    l,dl = state_vec
    v,ks,x,m = params
    return np.array([dl,((np.cosh(x/2.)**4. * v**2. * np.tanh(l))/np.cosh(l)**2. - (2. * ks * (l - x/2.))/m)])

# Setup system jacobian (Needed since method is implicit)
def dynjac_h3exactbar(state_vec, params):
    l,dl = state_vec
    v,ks,x,m = params
    return np.array([
        [0.,1.],
        [-(2. * ks/m) - v**2. * np.cosh(x/2)**4. * (-2. + np.cosh(2*l)) * 1./np.cosh(l)**4., 0.]
    ])

#####################
# H3 Sim Spring Bar #
#####################

def dynfunc_h3simbar(state_vec,params):
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
    v,ks,x,m = params

    # Distance Function
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # First derivatives of distance function
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)

    spa1 = (ks*(arccosh(d12) - x)*da1d12)/(m*sqrt(d12**2. - 1.))
    spb1 = (ks*(arccosh(d12) - x)*db1d12)/(m*sinh(a1)**2. * sqrt(d12**2. - 1.))
    spg1 = (ks*(arccosh(d12) - x)*dg1d12)/(m*sinh(a1)**2. * sin(b1)**2. * sqrt(d12**2. - 1.))
    spa2 = (ks*(arccosh(d12) - x)*da2d12)/(m*sqrt(d12**2. - 1.))
    spb2 = (ks*(arccosh(d12) - x)*db2d12)/(m*sinh(a2)**2. * sqrt(d12**2. - 1.))
    spg2 = (ks*(arccosh(d12) - x)*dg2d12)/(m*sinh(a2)**2. * sin(b2)**2. * sqrt(d12**2. - 1.))

    riga1 = 1./2.*sinh(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) - spa1
    rigb1 = -2. * ad1*bd1/tanh(a1) + .5*sin(2.*b1)*gd1**2 - spb1
    rigg1 = -2. * ad1*gd1/tanh(a1) - 2.*bd1*gd1/tan(b1) - spg1
    riga2 = 1./2.*sinh(2. * a2)*(bd2**2. + gd2**2. * sin(b2)**2.) - spa2
    rigb2 = -2. * ad2*bd2/tanh(a2) + .5*sin(2.*b2)*gd2**2 - spb2
    rigg2 = -2. * ad2*gd2/tanh(a2) - 2.*bd2*gd2/tan(b2) - spg2

    return np.array(
        [ad1,bd1,gd1,ad2,bd2,gd2,
        riga1,rigb1,rigg1,riga2,rigb2,rigg2])

# Setup system jacobian (Needed since method is implicit)
def dynjac_h3simbar(state_vec, params):
    def D12(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)
    
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

    def da1D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def db1D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*sin(b1)*sinh(a2)*cos(b2) - cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def dg1D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def da2D12a1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*sinh(a2) - cosh(a1)*cos(b1)*cosh(a2)*cos(b2) - cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

    def db2D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*cos(b1)*sinh(a2)*sin(b2) - cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

    def dg2D12a1(a1, b1, g1, a2, b2, g2):
        return -cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def db1D12b1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*cos(b1)*sinh(a2)*cos(b2) + sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def dg1D12b1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def db2D12b1(a1, b1, g1, a2, b2, g2):
        return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

    def dg2D12b1(a1, b1, g1, a2, b2, g2):
        return -sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def dg1D12g1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def dg2D12g1(a1, b1, g1, a2, b2, g2):
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

    def da2da1V12(m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1):
        # negative sign here so that the term can be added to geo terms later
        return -k/(m*f*sqrt( d12**2. - 1. ))*( (da1d12*da2d12)/sqrt( d12**2. - 1.) + ( arccosh(d12) - l )*( da2d12da1 - da1d12*(da2f/f + d12*da2d12/(d12**2. - 1.)) ) )

    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    v,ks,x,m = params

    # Distance function
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # First derivatives of distance function
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    # Second derivatives of distance function
    da1d12a1=da1D12a1(a1, b1, g1, a2, b2, g2)
    db1d12a1=db1D12a1(a1, b1, g1, a2, b2, g2)
    dg1d12a1=dg1D12a1(a1, b1, g1, a2, b2, g2)
    da2d12a1=da2D12a1(a1, b1, g1, a2, b2, g2)
    db2d12a1=db2D12a1(a1, b1, g1, a2, b2, g2)
    dg2d12a1=dg2D12a1(a1, b1, g1, a2, b2, g2)
    
    da1d12b1=db1d12a1
    db1d12b1=db1D12b1(a1, b1, g1, a2, b2, g2)
    dg1d12b1=dg1D12b1(a1, b1, g1, a2, b2, g2)
    da2d12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12b1=db2D12b1(a1, b1, g1, a2, b2, g2)
    dg2d12b1=dg2D12b1(a1, b1, g1, a2, b2, g2)

    da1d12g1=dg1d12a1
    db1d12g1=dg1d12b1
    dg1d12g1=dg1D12g1(a1, b1, g1, a2, b2, g2)
    da2d12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12g1=dg2D12g1(a1, b1, g1, a2, b2, g2)

    da1d12a2=da2d12a1
    db1d12a2=da2d12b1
    dg1d12a2=da2d12g1
    da2d12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
    db2d12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
    dg2d12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

    da1d12b2=db2d12a1
    db1d12b2=db2d12b1
    dg1d12b2=db2d12g1
    da2d12b2=db2d12a2
    db2d12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

    da1d12g2=dg2d12a1
    db1d12g2=dg2d12b1
    dg1d12g2=dg2d12g1
    da2d12g2=dg2d12a2
    db2d12g2=dg2d12b2
    dg2d12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)

    #---------- particle 1

    da1riga1 = cosh(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + da2da1V12(m, 1., ks, x, d12, da1d12, da1d12, 0., da1d12a1)
    db1riga1 = cos(b1)*sin(b1)*sinh(2.*a1)*gd1**2.           + da2da1V12(m, 1., ks, x, d12, da1d12, db1d12, 0., db1d12a1)
    dg1riga1 = 0.                                            + da2da1V12(m, 1., ks, x, d12, da1d12, dg1d12, 0., dg1d12a1)

    da2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, da2d12, 0., da2d12a1)
    db2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, db2d12, 0., db2d12a1)
    dg2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, dg2d12, 0., dg2d12a1)

    dad1riga1 = 0.                          + 0.
    dbd1riga1 = sinh(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sinh(2.*a1)*gd1 + 0.

    dad2riga1 = 0. + 0.
    dbd2riga1 = 0. + 0.
    dgd2riga1 = 0. + 0.

    #----------

    da1rigb1 = 2./sinh(a1)**2.*(ad1*bd1) + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, da1d12, sinh(2.*a1), da1d12b1)
    db1rigb1 = cos(2.*b1)*gd1**2.        + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, db1d12, 0.,          db1d12b1)
    dg1rigb1 = 0.                        + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, dg1d12, 0.,          dg1d12b1)

    da2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, da2d12, 0., da2d12b1)
    db2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, db2d12, 0., db2d12b1)
    dg2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, dg2d12, 0., dg2d12b1)

    dad1rigb1 = -2./tanh(a1)*bd1 + 0.
    dbd1rigb1 = -2./tanh(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1    + 0.

    dad2rigb1 = 0. + 0.
    dbd2rigb1 = 0. + 0.
    dgd2rigb1 = 0. + 0.

    #------------

    da1rigg1 = 2./sinh(a1)**2.*(ad1*gd1) + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da1d12, sinh(2.*a1)*sin(b1)**2., da1d12g1)
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)  + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db1d12, sin(2.*b1)*sinh(a1)**2., db1d12g1)
    dg1rigg1 = 0.                        + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg1d12, 0.,                      dg1d12g1)

    da2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da2d12, 0., da2d12g1)
    db2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db2d12, 0., db2d12g1)
    dg2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg2d12, 0., dg2d12g1)

    dad1rigg1 = -2./tanh(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                  + 0.
    dgd1rigg1 = -2.*(ad1/tanh(a1) + bd1/tan(b1)) + 0.

    dad2rigg1 = 0. + 0.
    dbd2rigg1 = 0. + 0.
    dgd2rigg1 = 0. + 0.

    #------------- particle 2

    da1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, da1d12, 0., da1d12a2)
    db1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, db1d12, 0., db1d12a2)
    dg1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, dg1d12, 0., dg1d12a2)

    da2riga2 = cosh(2.*a2) * (bd2**2. + sin(b2)**2.*gd2**2.) + da2da1V12(m, 1., ks, x, d12, da2d12, da2d12, 0., da2d12a2)
    db2riga2 = cos(b2)*sin(b2)*sinh(2.*a2)*gd2**2.           + da2da1V12(m, 1., ks, x, d12, da2d12, db2d12, 0., db2d12a2)
    dg2riga2 = 0.                                            + da2da1V12(m, 1., ks, x, d12, da2d12, dg2d12, 0., dg2d12a2)

    dad1riga2 = 0. + 0.
    dbd1riga2 = 0. + 0.
    dgd1riga2 = 0. + 0.

    dad2riga2 = 0.                          + 0.
    dbd2riga2 = sinh(2.*a2)*bd2             + 0.
    dgd2riga2 = sin(b2)**2.*sinh(2.*a2)*gd2 + 0.

    #----------

    da1rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, da1d12, 0., da1d12b2)
    db1rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, db1d12, 0., db1d12b2)
    dg1rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, dg1d12, 0., dg1d12b2)

    da2rigb2 = 2./sinh(a2)**2.*(ad2*bd2) + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, da2d12, sinh(2.*a2), da2d12b2)
    db2rigb2 = cos(2.*b2)*gd2**2.        + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, db2d12, 0.,          db2d12b2)
    dg2rigb2 = 0.                        + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, dg2d12, 0.,          dg2d12b2)

    dad1rigb2 = -2./tanh(a2)*bd2 + 0.
    dbd1rigb2 = -2./tanh(a2)*ad2 + 0.
    dgd1rigb2 = sin(2.*b2)*g2    + 0.

    dad2rigb2 = 0. + 0.
    dbd2rigb2 = 0. + 0.
    dgd2rigb2 = 0. + 0.

    #------------

    da1rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da1d12, 0., da1d12g2)
    db1rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db1d12, 0., db1d12g2)
    dg1rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg1d12, 0., dg1d12g2)

    da2rigg2 = 2./sinh(a2)**2.*(ad2*gd2) + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da2d12, sinh(2.*a2)*sin(b2)**2., da2d12g2)
    db2rigg2 = 2./sin(b2)**2.*(bd2*gd2)  + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db2d12, sin(2.*b2)*sinh(a2)**2., db2d12g2)
    dg2rigg2 = 0.                        + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg2d12, 0.,                      dg2d12g2)

    dad1rigg2 = -2./tanh(a2)*gd2                 + 0.
    dbd1rigg2 = -2./tan(b2)*gd2                  + 0.
    dgd1rigg2 = -2.*(ad2/tanh(a2) + bd2/tan(b2)) + 0.

    dad2rigg2 = 0. + 0.
    dbd2rigg2 = 0. + 0.
    dgd2rigg2 = 0. + 0.

    return np.array([
        [0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, da2riga1,db2riga1,dg2riga1, dad1riga1,dbd1riga1,dgd1riga1, dad2riga1,dbd2riga1,dgd2riga1],
        [da1rigb1,db1rigb1,dg1rigb1, da2rigb1,db2rigb1,dg2rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1, dad2rigb1,dbd2rigb1,dgd2rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, da2rigg1,db2rigg1,dg2rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1, dad2rigg1,dbd2rigg1,dgd2rigg1],
        [da1riga2,db1riga2,dg1riga2, da2riga2,db2riga2,dg2riga2, dad1riga2,dbd1riga2,dgd1riga2, dad2riga2,dbd2riga2,dgd2riga2],
        [da1rigb2,db1rigb2,dg1rigb2, da2rigb2,db2rigb2,dg2rigb2, dad1rigb2,dbd1rigb2,dgd1rigb2, dad2rigb2,dbd2rigb2,dgd2rigb2],
        [da1rigg2,db1rigg2,dg1rigg2, da2rigg2,db2rigg2,dg2rigg2, dad1rigg2,dbd1rigg2,dgd1rigg2, dad2rigg2,dbd2rigg2,dgd2rigg2]
    ])

############################
# H3 Exact Spring Triangle #
############################

def dynfunc_h3exacttriangle(state_vec,params):
    l,H,dl,dH = state_vec
    v,ks,x,m = params
    P = m*v*(1. + 2.*cosh(x/2.)**2.)
    aux1 = cosh(l)*cosh(H)
    return np.array([
        dl,dH,
        1./(2.*m**2.)*(-4.*ks*m*l + 2.*ks*m*(x + ((x - arccosh(aux1))*cosh(H)*sinh(l))/(sqrt(aux1**2. - 1.))) + (sinh(2.*l)*(P - m*dH)**2.)/(2. + cosh(2.*l))**2.),
        ((2.*ks*(x - arccosh(aux1))*cosh(l)*sinh(H))/(sqrt(aux1**2. - 1.)) + (2.*sinh(2.*l)*(P - m*dH)*dl)/(2. + cosh(2.*l))**2.)/(m - m/(2. + cosh(2.*l)))])

# Setup system jacobian (Needed since method is implicit)
def dynjac_h3exacttriangle(state_vec, params):
    l,H,dl,dH = state_vec
    v,ks,x,m = params
    P = m*v*(1. + 2.*cosh(x/2.)**2.)
    aux1 = cosh(l)*cosh(H)
    aux2 = sinh(l)*sinh(H)
    return np.array([
        [0.,0.,1.,0.],
        [0.,0.,0.,1.],
        [(1./(2.*m**2.))*(-4.*ks*m + (2.*(-dH*m + P)**2.*cosh(2.*l))/(2. + cosh(2.*l))**2. + 2.*ks*m*(((x - arccosh(aux1))*aux1)/(sqrt(aux1**2. - 1.)) - ((x - arccosh(aux1))*cosh(H)**2.*sinh(l)**2.)/(2.*sqrt(-1. + aux1)*(1. + aux1)**(3./2.)) - (cosh(H)**2.*sinh(l)**2.)/(-1 + aux1**2.) - ((x - arccosh(aux1))*cosh(H)**2.*sinh(l)**2.)/(2*(-1 + aux1)**(3./2.)*sqrt(1. + aux1))) - (4.*(-dH*m + P)**2.*sinh(2.*l)**2.)/(2. + cosh(2.*l))**3.),
         -((ks*(x - arccosh(aux1) + aux1*sqrt(-1 + aux1**2.))*aux2)/(m*(-1. + aux1**2.)**(3./2.))),
         0.,
         ((dH*m - P)*sinh(2.*l))/(m*(2. + cosh(2.*l))**2.)],
        [-((2.*m*sinh(2.*l)*((2.*ks*(x - arccosh(aux1))*cosh(l)*sinh(H))/(sqrt(-1. + aux1**2.)) + (2.*dl*(-dH*m + P)*sinh(2.*l))/(2. + cosh(2.*l))**2.))/((2 + cosh(2.*l))**2.*(m - m/(2. + cosh(2.*l)))**2.)) + (1./(m - m/(2. + cosh(2.*l))))*((4.*dl*(-dH*m + P)*cosh(2.*l))/(2. + cosh(2.*l))**2. - (ks*(x - arccosh(aux1))*aux1*aux2)/(sqrt(-1. + aux1)*(1. + aux1)**(3./2.)) - (2.*ks*aux1*aux2)/((-1 + aux1**2.)) - (ks*(x - arccosh(aux1))*aux1*aux2)/((-1. + aux1)**(3./2.)*sqrt(1. + aux1)) + (2.*ks*(x - arccosh(aux1))*aux2)/(sqrt(-1. + aux1**2.)) - (8.*dl*(-dH*m + P)*sinh(2.*l)**2.)/(2. + cosh(2.*l))**3.),
         -((ks*cosh(l)*(2. + cosh(2.*l))*(-x*cosh(H)**3. + sqrt(-1 + aux1**2.)/cosh(l)*sinh(H)**2. + x*cosh(H)*(1./cosh(l)**2. + sinh(H)**2.) + arccosh(aux1)*cosh(H)*tanh(l)**2.))/(m*(-1. + aux1**2.)**(3./2.))),
         (2.*(-dH*m + P)*tanh(l))/(m*(2. + cosh(2.*l))),
         -((2.*dl*tanh(l))/(2. + cosh(2.*l)))]
        ])

##########################
# H3 Sim Spring Triangle #
##########################

def dynfunc_h3simtriangle(state_vec,params):
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
    
    a1,b1,g1,a2,b2,g2,a3,b3,g3,ad1,bd1,gd1,ad2,bd2,gd2,ad3,bd3,gd3 = state_vec
    v,ks,x,m = params

    # Distance Functions
    # spring 12
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # spring 13
    d13 = D12(a1, b1, g1, a3, b3, g3)
    # spring 23
    d23 = D12(a2, b2, g2, a3, b3, g3)
    # First derivatives of distance function
    # spring 12
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    # spring 13
    da1d13 = da1D12(a1, b1, g1, a3, b3, g3)
    db1d13 = db1D12(a1, b1, g1, a3, b3, g3)
    dg1d13 = dg1D12(a1, b1, g1, a3, b3, g3)
    da3d13 = da1D12(a3, b3, g3, a1, b1, g1)
    db3d13 = db1D12(a3, b3, g3, a1, b1, g1)
    dg3d13 = dg1D12(a3, b3, g3, a1, b1, g1)
    # spring 23
    da2d23 = da1D12(a2, b2, g2, a3, b3, g3)
    db2d23 = db1D12(a2, b2, g2, a3, b3, g3)
    dg2d23 = dg1D12(a2, b2, g2, a3, b3, g3)
    da3d23 = da1D12(a3, b3, g3, a2, b2, g2)
    db3d23 = db1D12(a3, b3, g3, a2, b2, g2)
    dg3d23 = dg1D12(a3, b3, g3, a2, b2, g2)

    def da1V12(m, f, k, l, d12, da1d12):
        # negative sign here so that the term can be added to geo terms later
        return -(k*(arccosh(d12) - l)*da1d12)/(m*f*sqrt( d12**2. - 1. ))

    
    spa1 = da1V12(m,1.,ks,x,d12,da1d12)                         + da1V12(m,1.,ks,x,d13,da1d13)
    spb1 = da1V12(m,sinh(a1)**2.,ks,x,d12,db1d12)               + da1V12(m,sinh(a1)**2.,ks,x,d13,db1d13)
    spg1 = da1V12(m,sinh(a1)**2.*sin(b1)**2.,ks,x,d12,dg1d12)   + da1V12(m,sinh(a1)**2.*sin(b1)**2.,ks,x,d13,dg1d13)
    spa2 = da1V12(m,1.,ks,x,d12,da2d12)                         + da1V12(m,1.,ks,x,d23,da2d23)
    spb2 = da1V12(m,sinh(a2)**2.,ks,x,d12,db2d12)               + da1V12(m,sinh(a2)**2.,ks,x,d23,db2d23)
    spg2 = da1V12(m,sinh(a2)**2.*sin(b2)**2.,ks,x,d12,dg2d12)   + da1V12(m,sinh(a2)**2.*sin(b2)**2.,ks,x,d23,dg2d23)
    spa3 = da1V12(m,1.,ks,x,d13,da3d13)                         + da1V12(m,1.,ks,x,d23,da3d23)
    spb3 = da1V12(m,sinh(a3)**2.,ks,x,d13,db3d13)               + da1V12(m,sinh(a3)**2.,ks,x,d23,db3d23)
    spg3 = da1V12(m,sinh(a3)**2.*sin(b3)**2.,ks,x,d13,dg3d13)   + da1V12(m,sinh(a3)**2.*sin(b3)**2.,ks,x,d23,dg3d23)

    riga1 = 1./2.*sinh(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.)   + spa1
    rigb1 = -2. * ad1*bd1/tanh(a1) + .5*sin(2.*b1)*gd1**2           + spb1
    rigg1 = -2. * ad1*gd1/tanh(a1) - 2.*bd1*gd1/tan(b1)             + spg1
    riga2 = 1./2.*sinh(2. * a2)*(bd2**2. + gd2**2. * sin(b2)**2.)   + spa2
    rigb2 = -2. * ad2*bd2/tanh(a2) + .5*sin(2.*b2)*gd2**2           + spb2
    rigg2 = -2. * ad2*gd2/tanh(a2) - 2.*bd2*gd2/tan(b2)             + spg2
    riga3 = 1./2.*sinh(2. * a3)*(bd3**2. + gd3**2. * sin(b3)**2.)   + spa3
    rigb3 = -2. * ad3*bd3/tanh(a3) + .5*sin(2.*b3)*gd3**2           + spb3
    rigg3 = -2. * ad3*gd3/tanh(a3) - 2.*bd3*gd3/tan(b3)             + spg3

    return np.array(
        [ad1,bd1,gd1,ad2,bd2,gd2,ad3,bd3,gd3,
        riga1,rigb1,rigg1,riga2,rigb2,rigg2,riga3,rigb3,rigg3])

# Setup system jacobian (Needed since method is implicit)
def dynjac_h3simtriangle(state_vec, params):
    def D12(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)
    
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

    def da1D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*cosh(a2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2) - sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def db1D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*sin(b1)*sinh(a2)*cos(b2) - cosh(a1)*cos(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def dg1D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def da2D12a1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*sinh(a2) - cosh(a1)*cos(b1)*cosh(a2)*cos(b2) - cosh(a1)*sin(b1)*cosh(a2)*sin(b2)*cos(g1 - g2)

    def db2D12a1(a1, b1, g1, a2, b2, g2):
        return cosh(a1)*cos(b1)*sinh(a2)*sin(b2) - cosh(a1)*sin(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

    def dg2D12a1(a1, b1, g1, a2, b2, g2):
        return -cosh(a1)*sin(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def db1D12b1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*cos(b1)*sinh(a2)*cos(b2) + sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def dg1D12b1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def db2D12b1(a1, b1, g1, a2, b2, g2):
        return -sinh(a1)*sin(b1)*sinh(a2)*sin(b2) - sinh(a1)*cos(b1)*sinh(a2)*cos(b2)*cos(g1 - g2)

    def dg2D12b1(a1, b1, g1, a2, b2, g2):
        return -sinh(a1)*cos(b1)*sinh(a2)*sin(b2)*sin(g1 - g2)

    def dg1D12g1(a1, b1, g1, a2, b2, g2):
        return sinh(a1)*sin(b1)*sinh(a2)*sin(b2)*cos(g1 - g2)

    def dg2D12g1(a1, b1, g1, a2, b2, g2):
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

    def da2da1V12(m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1):
        # negative sign here so that the term can be added to geo terms later
        return -k/(m*f*sqrt( d12**2. - 1. ))*( (da1d12*da2d12)/sqrt( d12**2. - 1.) + ( arccosh(d12) - l )*( da2d12da1 - da1d12*(da2f/f + d12*da2d12/(d12**2. - 1.)) ) )

    a1,b1,g1,a2,b2,g2,a3,b3,g3,ad1,bd1,gd1,ad2,bd2,gd2,ad3,bd3,gd3 = state_vec
    v,ks,x,m = params

    # Distance Functions
    # spring 12
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # spring 13
    d13 = D12(a1, b1, g1, a3, b3, g3)
    # spring 23
    d23 = D12(a2, b2, g2, a3, b3, g3)
    # First derivatives of distance function
    # spring 12
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    # spring 13
    da1d13 = da1D12(a1, b1, g1, a3, b3, g3)
    db1d13 = db1D12(a1, b1, g1, a3, b3, g3)
    dg1d13 = dg1D12(a1, b1, g1, a3, b3, g3)
    da3d13 = da1D12(a3, b3, g3, a1, b1, g1)
    db3d13 = db1D12(a3, b3, g3, a1, b1, g1)
    dg3d13 = dg1D12(a3, b3, g3, a1, b1, g1)
    # spring 23
    da2d23 = da1D12(a2, b2, g2, a3, b3, g3)
    db2d23 = db1D12(a2, b2, g2, a3, b3, g3)
    dg2d23 = dg1D12(a2, b2, g2, a3, b3, g3)
    da3d23 = da1D12(a3, b3, g3, a2, b2, g2)
    db3d23 = db1D12(a3, b3, g3, a2, b2, g2)
    dg3d23 = dg1D12(a3, b3, g3, a2, b2, g2)
    # Second derivatives of distance function
    # spring 12
    da1d12a1=da1D12a1(a1, b1, g1, a2, b2, g2)
    db1d12a1=db1D12a1(a1, b1, g1, a2, b2, g2)
    dg1d12a1=dg1D12a1(a1, b1, g1, a2, b2, g2)
    da2d12a1=da2D12a1(a1, b1, g1, a2, b2, g2)
    db2d12a1=db2D12a1(a1, b1, g1, a2, b2, g2)
    dg2d12a1=dg2D12a1(a1, b1, g1, a2, b2, g2)
    
    da1d12b1=db1d12a1
    db1d12b1=db1D12b1(a1, b1, g1, a2, b2, g2)
    dg1d12b1=dg1D12b1(a1, b1, g1, a2, b2, g2)
    da2d12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12b1=db2D12b1(a1, b1, g1, a2, b2, g2)
    dg2d12b1=dg2D12b1(a1, b1, g1, a2, b2, g2)

    da1d12g1=dg1d12a1
    db1d12g1=dg1d12b1
    dg1d12g1=dg1D12g1(a1, b1, g1, a2, b2, g2)
    da2d12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12g1=dg2D12g1(a1, b1, g1, a2, b2, g2)

    da1d12a2=da2d12a1
    db1d12a2=da2d12b1
    dg1d12a2=da2d12g1
    da2d12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
    db2d12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
    dg2d12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

    da1d12b2=db2d12a1
    db1d12b2=db2d12b1
    dg1d12b2=db2d12g1
    da2d12b2=db2d12a2
    db2d12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

    da1d12g2=dg2d12a1
    db1d12g2=dg2d12b1
    dg1d12g2=dg2d12g1
    da2d12g2=dg2d12a2
    db2d12g2=dg2d12b2
    dg2d12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)
    # spring 13
    da1d13a1=da1D12a1(a1, b1, g1, a3, b3, g3)
    db1d13a1=db1D12a1(a1, b1, g1, a3, b3, g3)
    dg1d13a1=dg1D12a1(a1, b1, g1, a3, b3, g3)
    da3d13a1=da2D12a1(a1, b1, g1, a3, b3, g3)
    db3d13a1=db2D12a1(a1, b1, g1, a3, b3, g3)
    dg3d13a1=dg2D12a1(a1, b1, g1, a3, b3, g3)
    
    da1d13b1=db1d13a1
    db1d13b1=db1D12b1(a1, b1, g1, a3, b3, g3)
    dg1d13b1=dg1D12b1(a1, b1, g1, a3, b3, g3)
    da3d13b1 = db2D12a1(a3, b3, g3, a1, b1, g1)
    db3d13b1=db2D12b1(a1, b1, g1, a3, b3, g3)
    dg3d13b1=dg2D12b1(a1, b1, g1, a3, b3, g3)

    da1d13g1=dg1d13a1
    db1d13g1=dg1d13b1
    dg1d13g1=dg1D12g1(a1, b1, g1, a3, b3, g3)
    da3d13g1 = dg2D12a1(a3, b3, g3, a1, b1, g1)
    db3d13g1 = dg2D12b1(a3, b3, g3, a1, b1, g1)
    dg3d13g1=dg2D12g1(a1, b1, g1, a3, b3, g3)

    da1d13a3=da3d13a1
    db1d13a3=da3d13b1
    dg1d13a3=da3d13g1
    da3d13a3 = da1D12a1(a3, b3, g3, a1, b1, g1)
    db3d13a3 = db1D12a1(a3, b3, g3, a1, b1, g1)
    dg3d13a3 = dg1D12a1(a3, b3, g3, a1, b1, g1)

    da1d13b3=db3d13a1
    db1d13b3=db3d13b1
    dg1d13b3=db3d13g1
    da3d13b3=db3d13a3
    db3d13b3 = db1D12b1(a3, b3, g3, a1, b1, g1)
    dg3d13b3 = dg1D12b1(a3, b3, g3, a1, b1, g1)

    da1d13g3=dg3d13a1
    db1d13g3=dg3d13b1
    dg1d13g3=dg3d13g1
    da3d13g3=dg3d13a3
    db3d13g3=dg3d13b3
    dg3d13g3 = dg1D12g1(a3, b3, g3, a1, b1, g1)
    # spring 23
    da2d23a2=da1D12a1(a2, b2, g2, a3, b3, g3)
    db2d23a2=db1D12a1(a2, b2, g2, a3, b3, g3)
    dg2d23a2=dg1D12a1(a2, b2, g2, a3, b3, g3)
    da3d23a2=da2D12a1(a2, b2, g2, a3, b3, g3)
    db3d23a2=db2D12a1(a2, b2, g2, a3, b3, g3)
    dg3d23a2=dg2D12a1(a2, b2, g2, a3, b3, g3)

    da2d23b2=db2d23a2
    db2d23b2=db1D12b1(a2, b2, g2, a3, b3, g3)
    dg2d23b2=dg1D12b1(a2, b2, g2, a3, b3, g3)
    da3d23b2 = db2D12a1(a3, b3, g3, a2, b2, g2)
    db3d23b2=db2D12b1(a2, b2, g2, a3, b3, g3)
    dg3d23b2=dg2D12b1(a2, b2, g2, a3, b3, g3)

    da2d23g2=dg2d23a2
    db2d23g2=dg2d23b2
    dg2d23g2=dg1D12g1(a2, b2, g2, a3, b3, g3)
    da3d23g2 = dg2D12a1(a3, b3, g3, a2, b2, g2)
    db3d23g2 = dg2D12b1(a3, b3, g3, a2, b2, g2)
    dg3d23g2=dg2D12g1(a2, b2, g2, a3, b3, g3)

    da2d23a3=da3d23a2
    db2d23a3=da3d23b2
    dg2d23a3=da3d23g2
    da3d23a3 = da1D12a1(a3, b3, g3, a2, b2, g2)
    db3d23a3 = db1D12a1(a3, b3, g3, a2, b2, g2)
    dg3d23a3 = dg1D12a1(a3, b3, g3, a2, b2, g2)

    da2d23b3=db3d23a2
    db2d23b3=db3d23b2
    dg2d23b3=db3d23g2
    da3d23b3=db3d23a3
    db3d23b3 = db1D12b1(a3, b3, g3, a2, b2, g2)
    dg3d23b3 = dg1D12b1(a3, b3, g3, a2, b2, g2)

    da2d23g3=dg3d23a2
    db2d23g3=dg3d23b2
    dg2d23g3=dg3d23g2
    da3d23g3=dg3d23a3
    db3d23g3=dg3d23b3
    dg3d23g3 = dg1D12g1(a3, b3, g3, a2, b2, g2)

    #---------- particle 1

    da1riga1 = cosh(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + da2da1V12(m, 1., ks, x, d12, da1d12, da1d12, 0., da1d12a1) + da2da1V12(m, 1., ks, x, d13, da1d13, da1d13, 0., da1d13a1)
    db1riga1 = cos(b1)*sin(b1)*sinh(2.*a1)*gd1**2.           + da2da1V12(m, 1., ks, x, d12, da1d12, db1d12, 0., db1d12a1) + da2da1V12(m, 1., ks, x, d13, da1d13, db1d13, 0., db1d13a1)
    dg1riga1 = 0.                                            + da2da1V12(m, 1., ks, x, d12, da1d12, dg1d12, 0., dg1d12a1) + da2da1V12(m, 1., ks, x, d13, da1d13, dg1d13, 0., dg1d13a1)

    da2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, da2d12, 0., da2d12a1)
    db2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, db2d12, 0., db2d12a1)
    dg2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, dg2d12, 0., dg2d12a1)

    da3riga1 = 0. + da2da1V12(m, 1., ks, x, d13, da1d12, da3d13, 0., da3d13a1)
    db3riga1 = 0. + da2da1V12(m, 1., ks, x, d13, da1d12, db3d13, 0., db3d13a1)
    dg3riga1 = 0. + da2da1V12(m, 1., ks, x, d13, da1d12, dg3d13, 0., dg3d13a1)

    dad1riga1 = 0.                          + 0.
    dbd1riga1 = sinh(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sinh(2.*a1)*gd1 + 0.

    dad2riga1 = 0. + 0.
    dbd2riga1 = 0. + 0.
    dgd2riga1 = 0. + 0.

    dad3riga1 = 0. + 0.
    dbd3riga1 = 0. + 0.
    dgd3riga1 = 0. + 0.

    #----------

    da1rigb1 = 2./sinh(a1)**2.*(ad1*bd1) + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, da1d12, sinh(2.*a1), da1d12b1) + da2da1V12(m, sinh(a1)**2., ks, x, d13, db1d13, da1d13, sinh(2.*a1), da1d13b1)
    db1rigb1 = cos(2.*b1)*gd1**2.        + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, db1d12, 0.,          db1d12b1) + da2da1V12(m, sinh(a1)**2., ks, x, d13, db1d13, db1d13, 0.,          db1d13b1)
    dg1rigb1 = 0.                        + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, dg1d12, 0.,          dg1d12b1) + da2da1V12(m, sinh(a1)**2., ks, x, d13, db1d13, dg1d13, 0.,          dg1d13b1)

    da2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, da2d12, 0., da2d12b1)
    db2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, db2d12, 0., db2d12b1)
    dg2rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d12, db1d12, dg2d12, 0., dg2d12b1)

    da3rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d13, db1d13, da3d13, 0., da3d13b1)
    db3rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d13, db1d13, db3d13, 0., db3d13b1)
    dg3rigb1 = 0. + da2da1V12(m, sinh(a1)**2., ks, x, d13, db1d13, dg3d13, 0., dg3d13b1)

    dad1rigb1 = -2./tanh(a1)*bd1 + 0.
    dbd1rigb1 = -2./tanh(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1    + 0.

    dad2rigb1 = 0. + 0.
    dbd2rigb1 = 0. + 0.
    dgd2rigb1 = 0. + 0.

    dad3rigb1 = 0. + 0.
    dbd3rigb1 = 0. + 0.
    dgd3rigb1 = 0. + 0.

    #------------

    da1rigg1 = 2./sinh(a1)**2.*(ad1*gd1) + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da1d12, sinh(2.*a1)*sin(b1)**2., da1d12g1) + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, da1d13, sinh(2.*a1)*sin(b1)**2., da1d13g1)
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)  + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db1d12, sin(2.*b1)*sinh(a1)**2., db1d12g1) + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, db1d13, sin(2.*b1)*sinh(a1)**2., db1d13g1)
    dg1rigg1 = 0.                        + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg1d12, 0.,                      dg1d12g1) + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, dg1d13, 0.,                      dg1d13g1)

    da2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da2d12, 0., da2d12g1)
    db2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db2d12, 0., db2d12g1)
    dg2rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg2d12, 0., dg2d12g1)

    da3rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, da3d13, 0., da3d13g1)
    db3rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, db3d13, 0., db3d13g1)
    dg3rigg1 = 0. + da2da1V12(m, sinh(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, dg3d13, 0., dg3d13g1)

    dad1rigg1 = -2./tanh(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                  + 0.
    dgd1rigg1 = -2.*(ad1/tanh(a1) + bd1/tan(b1)) + 0.

    dad2rigg1 = 0. + 0.
    dbd2rigg1 = 0. + 0.
    dgd2rigg1 = 0. + 0.

    dad3rigg1 = 0. + 0.
    dbd3rigg1 = 0. + 0.
    dgd3rigg1 = 0. + 0.

    #------------- particle 2

    da1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, da1d12, 0., da1d12a2)
    db1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, db1d12, 0., db1d12a2)
    dg1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, dg1d12, 0., dg1d12a2)

    da2riga2 = cosh(2.*a2) * (bd2**2. + sin(b2)**2.*gd2**2.) + da2da1V12(m, 1., ks, x, d12, da2d12, da2d12, 0., da2d12a2) + da2da1V12(m, 1., ks, x, d23, da2d23, da2d23, 0., da2d23a2)
    db2riga2 = cos(b2)*sin(b2)*sinh(2.*a2)*gd2**2.           + da2da1V12(m, 1., ks, x, d12, da2d12, db2d12, 0., db2d12a2) + da2da1V12(m, 1., ks, x, d23, da2d23, db2d23, 0., db2d23a2)
    dg2riga2 = 0.                                            + da2da1V12(m, 1., ks, x, d12, da2d12, dg2d12, 0., dg2d12a2) + da2da1V12(m, 1., ks, x, d23, da2d23, dg2d23, 0., dg2d23a2)

    da3riga2 = 0. + da2da1V12(m, 1., ks, x, d23, da2d23, da3d23, 0., da3d23a2)
    db3riga2 = 0. + da2da1V12(m, 1., ks, x, d23, da2d23, db3d23, 0., db3d23a2)
    dg3riga2 = 0. + da2da1V12(m, 1., ks, x, d23, da2d23, dg3d23, 0., dg3d23a2)
    
    dad1riga2 = 0. + 0.
    dbd1riga2 = 0. + 0.
    dgd1riga2 = 0. + 0.

    dad2riga2 = 0.                          + 0.
    dbd2riga2 = sinh(2.*a2)*bd2             + 0.
    dgd2riga2 = sin(b2)**2.*sinh(2.*a2)*gd2 + 0.

    dad3riga2 = 0. + 0.
    dbd3riga2 = 0. + 0.
    dgd3riga2 = 0. + 0.

    #----------

    da1rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, da1d12, 0., da1d12b2)
    db1rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, db1d12, 0., db1d12b2)
    dg1rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, dg1d12, 0., dg1d12b2)

    da2rigb2 = 2./sinh(a2)**2.*(ad2*bd2) + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, da2d12, sinh(2.*a2), da2d12b2) + da2da1V12(m, sinh(a2)**2., ks, x, d23, db2d23, da2d23, sinh(2.*a2), da2d23b2)
    db2rigb2 = cos(2.*b2)*gd2**2.        + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, db2d12, 0.,          db2d12b2) + da2da1V12(m, sinh(a2)**2., ks, x, d23, db2d23, db2d23, 0.,          db2d23b2)
    dg2rigb2 = 0.                        + da2da1V12(m, sinh(a2)**2., ks, x, d12, db2d12, dg2d12, 0.,          dg2d12b2) + da2da1V12(m, sinh(a2)**2., ks, x, d23, db2d23, dg2d23, 0.,          dg2d23b2)

    da3rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d23, db2d23, da3d23, 0., da3d23b2)
    db3rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d23, db2d23, db3d23, 0., db3d23b2)
    dg3rigb2 = 0. + da2da1V12(m, sinh(a2)**2., ks, x, d23, db2d23, dg3d23, 0., dg3d23b2)
    
    dad1rigb2 = -2./tanh(a2)*bd2 + 0.
    dbd1rigb2 = -2./tanh(a2)*ad2 + 0.
    dgd1rigb2 = sin(2.*b2)*g2    + 0.

    dad2rigb2 = 0. + 0.
    dbd2rigb2 = 0. + 0.
    dgd2rigb2 = 0. + 0.

    dad3rigb2 = 0. + 0.
    dbd3rigb2 = 0. + 0.
    dgd3rigb2 = 0. + 0.

    #------------

    da1rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da1d12, 0., da1d12g2)
    db1rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db1d12, 0., db1d12g2)
    dg1rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg1d12, 0., dg1d12g2)

    da2rigg2 = 2./sinh(a2)**2.*(ad2*gd2) + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da2d12, sinh(2.*a2)*sin(b2)**2., da2d12g2) + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, da2d23, sinh(2.*a2)*sin(b2)**2., da2d23g2)
    db2rigg2 = 2./sin(b2)**2.*(bd2*gd2)  + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db2d12, sin(2.*b2)*sinh(a2)**2., db2d12g2) + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, db2d23, sin(2.*b2)*sinh(a2)**2., db2d23g2)
    dg2rigg2 = 0.                        + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg2d12, 0.,                      dg2d12g2) + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, dg2d23, 0.,                      dg2d23g2)

    da3rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, da3d23, 0., da3d23g2)
    db3rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, db3d23, 0., db3d23g2)
    dg3rigg2 = 0. + da2da1V12(m, sinh(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, dg3d23, 0., dg3d23g2)
    
    dad1rigg2 = -2./tanh(a2)*gd2                 + 0.
    dbd1rigg2 = -2./tan(b2)*gd2                  + 0.
    dgd1rigg2 = -2.*(ad2/tanh(a2) + bd2/tan(b2)) + 0.

    dad2rigg2 = 0. + 0.
    dbd2rigg2 = 0. + 0.
    dgd2rigg2 = 0. + 0.

    dad3rigg2 = 0. + 0.
    dbd3rigg2 = 0. + 0.
    dgd3rigg2 = 0. + 0.

    #------------- particle 3

    da1riga3 = 0. + da2da1V12(m, 1., ks, x, d13, da3d13, da1d13, 0., da1d13a3)
    db1riga3 = 0. + da2da1V12(m, 1., ks, x, d13, da3d13, db1d13, 0., db1d13a3)
    dg1riga3 = 0. + da2da1V12(m, 1., ks, x, d13, da3d13, dg1d13, 0., dg1d13a3)

    da2riga3 = 0. + da2da1V12(m, 1., ks, x, d23, da3d23, da2d23, 0., da2d23a3)
    db2riga3 = 0. + da2da1V12(m, 1., ks, x, d23, da3d23, db2d23, 0., db2d23a3)
    dg2riga3 = 0. + da2da1V12(m, 1., ks, x, d23, da3d23, dg2d23, 0., dg2d23a3)

    da3riga3 = cosh(2.*a3) * (bd3**2. + sin(b3)**2.*gd3**2.) + da2da1V12(m, 1., ks, x, d12, da3d13, da3d13, 0., da3d13a3) + da2da1V12(m, 1., ks, x, d23, da3d23, da3d23, 0., da3d23a3)
    db3riga3 = cos(b3)*sin(b3)*sinh(2.*a3)*gd3**2.           + da2da1V12(m, 1., ks, x, d12, da3d13, db3d13, 0., db3d13a3) + da2da1V12(m, 1., ks, x, d23, da3d23, db3d23, 0., db3d23a3)
    dg3riga3 = 0.                                            + da2da1V12(m, 1., ks, x, d12, da3d13, dg3d13, 0., dg3d13a3) + da2da1V12(m, 1., ks, x, d23, da3d23, dg3d23, 0., dg3d23a3)
    
    dad1riga3 = 0. + 0.
    dbd1riga3 = 0. + 0.
    dgd1riga3 = 0. + 0.

    dad2riga3 = 0. + 0.
    dbd2riga3 = 0. + 0.
    dgd2riga3 = 0. + 0.

    dad3riga3 = 0.                          + 0.
    dbd3riga3 = sinh(2.*a3)*bd3             + 0.
    dgd3riga3 = sin(b3)**2.*sinh(2.*a3)*gd3 + 0.

    #----------

    da1rigb3 = 0. + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, da1d13, 0., da1d13b3)
    db1rigb3 = 0. + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, db1d13, 0., db1d13b3)
    dg1rigb3 = 0. + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, dg1d13, 0., dg1d13b3)

    da2rigb3 = 0. + da2da1V12(m, sinh(a3)**2., ks, x, d23, db3d23, da2d23, 0., da2d23b3)
    db2rigb3 = 0. + da2da1V12(m, sinh(a3)**2., ks, x, d23, db3d23, db2d23, 0., db2d23b3)
    dg2rigb3 = 0. + da2da1V12(m, sinh(a3)**2., ks, x, d23, db3d23, dg2d23, 0., dg2d23b3)

    da3rigb3 = 2./sinh(a3)**2.*(ad3*bd3) + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, da3d13, sinh(2.*a3), da3d13b3) + da2da1V12(m, sinh(a3)**2., ks, x, d23, db3d23, da3d23, sinh(2.*a3), da3d23b3)
    db3rigb3 = cos(2.*b3)*gd3**2.        + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, db3d13, 0.,          db3d13b3) + da2da1V12(m, sinh(a3)**2., ks, x, d23, db3d23, db3d23, 0.,          db3d23b3)
    dg3rigb3 = 0.                        + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, dg3d13, 0.,          dg3d13b3) + da2da1V12(m, sinh(a3)**2., ks, x, d23, db3d23, dg3d23, 0.,          dg3d23b3)

    dad1rigb3 = 0. + 0.
    dbd1rigb3 = 0. + 0.
    dgd1rigb3 = 0. + 0.

    dad2rigb3 = 0. + 0.
    dbd2rigb3 = 0. + 0.
    dgd2rigb3 = 0. + 0.

    dad3rigb3 = -2./tanh(a3)*bd3 + 0.
    dbd3rigb3 = -2./tanh(a3)*ad3 + 0.
    dgd3rigb3 = sin(2.*b3)*g3    + 0.

    #------------

    da1rigg3 = 0. + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, da1d13, 0., da1d13g3)
    db1rigg3 = 0. + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, db1d13, 0., db1d13g3)
    dg1rigg3 = 0. + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, dg1d13, 0., dg1d13g3)

    da2rigg3 = 0. + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, da2d23, 0., da2d23g3)
    db2rigg3 = 0. + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, db2d23, 0., db2d23g3)
    dg2rigg3 = 0. + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, dg2d23, 0., dg2d23g3)

    da3rigg3 = 2./sinh(a3)**2.*(ad3*gd3) + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, da3d13, sinh(2.*a3)*sin(b3)**2., da3d13g3) + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, da3d23, sinh(2.*a3)*sin(b3)**2., da3d23g3)
    db3rigg3 = 2./sin(b3)**2.*(bd3*gd3)  + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, db3d13, sin(2.*b3)*sinh(a3)**2., db3d13g3) + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, db3d23, sin(2.*b3)*sinh(a3)**2., db3d23g3)
    dg3rigg3 = 0.                        + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, dg3d13, 0.,                      dg3d13g3) + da2da1V12(m, sinh(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, dg3d23, 0.,                      dg3d23g3)

    dad1rigg3 = 0. + 0.
    dbd1rigg3 = 0. + 0.
    dgd1rigg3 = 0. + 0.

    dad2rigg3 = 0. + 0.
    dbd2rigg3 = 0. + 0.
    dgd2rigg3 = 0. + 0.

    dad3rigg3 = -2./tanh(a3)*gd3                 + 0.
    dbd3rigg3 = -2./tan(b3)*gd3                  + 0.
    dgd3rigg3 = -2.*(ad3/tanh(a3) + bd3/tan(b3)) + 0.

    return np.array([
        [0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.],
        [da1riga1,db1riga1,dg1riga1, da2riga1,db2riga1,dg2riga1, da3riga1,db3riga1,dg3riga1, dad1riga1,dbd1riga1,dgd1riga1, dad2riga1,dbd2riga1,dgd2riga1, dad3riga1,dbd3riga1,dgd3riga1],
        [da1rigb1,db1rigb1,dg1rigb1, da2rigb1,db2rigb1,dg2rigb1, da3rigb1,db3rigb1,dg3rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1, dad2rigb1,dbd2rigb1,dgd2rigb1, dad3rigb1,dbd3rigb1,dgd3rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, da2rigg1,db2rigg1,dg2rigg1, da3rigg1,db3rigg1,dg3rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1, dad2rigg1,dbd2rigg1,dgd2rigg1, dad3rigg1,dbd3rigg1,dgd3rigg1],
        [da1riga2,db1riga2,dg1riga2, da2riga2,db2riga2,dg2riga2, da3riga2,db3riga2,dg3riga2, dad1riga2,dbd1riga2,dgd1riga2, dad2riga2,dbd2riga2,dgd2riga2, dad3riga2,dbd3riga2,dgd3riga2],
        [da1rigb2,db1rigb2,dg1rigb2, da2rigb2,db2rigb2,dg2rigb2, da3rigb2,db3rigb2,dg3rigb2, dad1rigb2,dbd1rigb2,dgd1rigb2, dad2rigb2,dbd2rigb2,dgd2rigb2, dad3rigb2,dbd3rigb2,dgd3rigb2],
        [da1rigg2,db1rigg2,dg1rigg2, da2rigg2,db2rigg2,dg2rigg2, da3rigg2,db3rigg2,dg3rigg2, dad1rigg2,dbd1rigg2,dgd1rigg2, dad2rigg2,dbd2rigg2,dgd2rigg2, dad3rigg2,dbd3rigg2,dgd3rigg2],
        [da1riga3,db1riga3,dg1riga3, da2riga3,db2riga3,dg2riga3, da3riga3,db3riga3,dg3riga3, dad1riga3,dbd1riga3,dgd1riga3, dad2riga3,dbd2riga3,dgd2riga3, dad3riga3,dbd3riga3,dgd3riga3],
        [da1rigb3,db1rigb3,dg1rigb3, da2rigb3,db2rigb3,dg2rigb3, da3rigb3,db3rigb3,dg3rigb3, dad1rigb3,dbd1rigb3,dgd1rigb3, dad2rigb3,dbd2rigb3,dgd2rigb3, dad3rigb3,dbd3rigb3,dgd3rigb3],
        [da1rigg3,db1rigg3,dg1rigg3, da2rigg3,db2rigg3,dg2rigg3, da3rigg3,db3rigg3,dg3rigg3, dad1rigg3,dbd1rigg3,dgd1rigg3, dad2rigg3,dbd2rigg3,dgd2rigg3, dad3rigg3,dbd3rigg3,dgd3rigg3]
    ])

############################
# H3 Sim Central Potential #
############################

def dynfunc_h3simcenpot(state_vec,params):
    
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,G,ms,m = params

    # Central Potential first derivative
    def h3da1cenpot(state_vec,params):
        a1,b1,g1,ad1,bd1,gd1 = state_vec
        v,G,ms,m = params
        return -G*ms/sinh(a1)**2.

    riga1 = 1./2.*sinh(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) + h3da1cenpot(state_vec,params)
    rigb1 = -2. * ad1*bd1/tanh(a1) + .5*sin(2.*b1)*gd1**2         + 0.
    rigg1 = -2. * ad1*gd1/tanh(a1) - 2.*bd1*gd1/tan(b1)           + 0.

    return np.array(
        [ad1,bd1,gd1,
        riga1,rigb1,rigg1])

# Setup system jacobian (Needed since method is implicit)
def dynjac_h3simcenpot(state_vec, params):
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,G,ms,m = params

    # Central Potential Second derivative
    def h3da1cenpotda1(state_vec,params):
        a1,b1,g1,ad1,bd1,gd1 = state_vec
        v,G,ms,m = params
        return 2.*G*ms/sinh(a1)**2./tanh(a1)**2.

    #---------- particle 1

    da1riga1 = cosh(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + h3da1cenpotda1(state_vec,params)
    db1riga1 = cos(b1)*sin(b1)*sinh(2.*a1)*gd1**2.           + 0.
    dg1riga1 = 0.                                            + 0.

    dad1riga1 = 0.                          + 0.
    dbd1riga1 = sinh(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sinh(2.*a1)*gd1 + 0.

    #----------

    da1rigb1 = 2./sinh(a1)**2.*(ad1*bd1) + 0.
    db1rigb1 = cos(2.*b1)*gd1**2.        + 0.
    dg1rigb1 = 0.                        + 0.

    dad1rigb1 = -2./tanh(a1)*bd1 + 0.
    dbd1rigb1 = -2./tanh(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1    + 0.

    #------------

    da1rigg1 = 2./sinh(a1)**2.*(ad1*gd1) + 0.
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)  + 0.
    dg1rigg1 = 0.                        + 0.

    dad1rigg1 = -2./tanh(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                  + 0.
    dgd1rigg1 = -2.*(ad1/tanh(a1) + bd1/tan(b1)) + 0.


    return np.array([
        [0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, dad1riga1,dbd1riga1,dgd1riga1],
        [da1rigb1,db1rigb1,dg1rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1]
    ])

#############################
# H3 Sim 2 Sphere Collision #
#############################

def dynfunc_h3sim2ballcol(state_vec,params):
    
    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    m1,m2,r1,r2 = params

    riga1 = 1./2.*sinh(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.)
    rigb1 = -2. * ad1*bd1/tanh(a1) + .5*sin(2.*b1)*gd1**2
    rigg1 = -2. * ad1*gd1/tanh(a1) - 2.*bd1*gd1/tan(b1)

    riga2 = 1./2.*sinh(2. * a2)*(bd2**2. + gd2**2. * sin(b2)**2.)
    rigb2 = -2. * ad2*bd2/tanh(a2) + .5*sin(2.*b2)*gd2**2
    rigg2 = -2. * ad2*gd2/tanh(a2) - 2.*bd2*gd2/tan(b2)

    return np.array(
        [ad1,bd1,gd1,ad2,bd2,gd2,
        riga1,rigb1,rigg1,riga2,rigb2,rigg2])

# Setup system jacobian (Needed since method is implicit)
def dynjac_h3sim2ballcol(state_vec, params):
    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    m1,m2,r1,r2 = params

    #---------- particle 1

    da1riga1 = cosh(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.)
    db1riga1 = cos(b1)*sin(b1)*sinh(2.*a1)*gd1**2.          
    dg1riga1 = 0.                                           

    da2riga1 = 0.
    db2riga1 = 0.
    dg2riga1 = 0.

    dad1riga1 = 0.                         
    dbd1riga1 = sinh(2.*a1)*bd1            
    dgd1riga1 = sin(b1)**2.*sinh(2.*a1)*gd1

    dad2riga1 = 0.
    dbd2riga1 = 0.
    dgd2riga1 = 0.

    #----------

    da1rigb1 = 2./sinh(a1)**2.*(ad1*bd1)
    db1rigb1 = cos(2.*b1)*gd1**2.       
    dg1rigb1 = 0.                       

    da2rigb1 = 0.
    db2rigb1 = 0.
    dg2rigb1 = 0.

    dad1rigb1 = -2./tanh(a1)*bd1
    dbd1rigb1 = -2./tanh(a1)*ad1
    dgd1rigb1 = sin(2.*b1)*g1   

    dad2rigb1 = 0.
    dbd2rigb1 = 0.
    dgd2rigb1 = 0.

    #------------

    da1rigg1 = 2./sinh(a1)**2.*(ad1*gd1)
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1) 
    dg1rigg1 = 0.                       

    da2rigg1 = 0.
    db2rigg1 = 0.
    dg2rigg1 = 0.

    dad1rigg1 = -2./tanh(a1)*gd1                
    dbd1rigg1 = -2./tan(b1)*gd1                 
    dgd1rigg1 = -2.*(ad1/tanh(a1) + bd1/tan(b1))

    dad2rigg1 = 0.
    dbd2rigg1 = 0.
    dgd2rigg1 = 0.

    #------------- particle 2

    da1riga2 = 0.
    db1riga2 = 0.
    dg1riga2 = 0.

    da2riga2 = cosh(2.*a2) * (bd2**2. + sin(b2)**2.*gd2**2.)
    db2riga2 = cos(b2)*sin(b2)*sinh(2.*a2)*gd2**2.          
    dg2riga2 = 0.                                           

    dad1riga2 = 0.
    dbd1riga2 = 0.
    dgd1riga2 = 0.

    dad2riga2 = 0.                         
    dbd2riga2 = sinh(2.*a2)*bd2            
    dgd2riga2 = sin(b2)**2.*sinh(2.*a2)*gd2

    #----------

    da1rigb2 = 0.
    db1rigb2 = 0.
    dg1rigb2 = 0.

    da2rigb2 = 2./sinh(a2)**2.*(ad2*bd2)
    db2rigb2 = cos(2.*b2)*gd2**2.       
    dg2rigb2 = 0.                       

    dad1rigb2 = -2./tanh(a2)*bd2
    dbd1rigb2 = -2./tanh(a2)*ad2
    dgd1rigb2 = sin(2.*b2)*g2   

    dad2rigb2 = 0.
    dbd2rigb2 = 0.
    dgd2rigb2 = 0.

    #------------

    da1rigg2 = 0.
    db1rigg2 = 0.
    dg1rigg2 = 0.

    da2rigg2 = 2./sinh(a2)**2.*(ad2*gd2)
    db2rigg2 = 2./sin(b2)**2.*(bd2*gd2) 
    dg2rigg2 = 0.                       

    dad1rigg2 = -2./tanh(a2)*gd2                
    dbd1rigg2 = -2./tan(b2)*gd2                 
    dgd1rigg2 = -2.*(ad2/tanh(a2) + bd2/tan(b2))

    dad2rigg2 = 0.
    dbd2rigg2 = 0.
    dgd2rigg2 = 0.

    return np.array([
        [0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, da2riga1,db2riga1,dg2riga1, dad1riga1,dbd1riga1,dgd1riga1, dad2riga1,dbd2riga1,dgd2riga1],
        [da1rigb1,db1rigb1,dg1rigb1, da2rigb1,db2rigb1,dg2rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1, dad2rigb1,dbd2rigb1,dgd2rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, da2rigg1,db2rigg1,dg2rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1, dad2rigg1,dbd2rigg1,dgd2rigg1],
        [da1riga2,db1riga2,dg1riga2, da2riga2,db2riga2,dg2riga2, dad1riga2,dbd1riga2,dgd1riga2, dad2riga2,dbd2riga2,dgd2riga2],
        [da1rigb2,db1rigb2,dg1rigb2, da2rigb2,db2rigb2,dg2rigb2, dad1rigb2,dbd1rigb2,dgd1rigb2, dad2rigb2,dbd2rigb2,dgd2rigb2],
        [da1rigg2,db1rigg2,dg1rigg2, da2rigg2,db2rigg2,dg2rigg2, dad1rigg2,dbd1rigg2,dgd1rigg2, dad2rigg2,dbd2rigg2,dgd2rigg2]
    ])





#########################################################
#################### S3 Test Systems ####################
#########################################################

###################
# S3 Sim Geodesic #
###################

def dynfunc_s3simgeo(state_vec,params):
    
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,m = params

    riga1 = 1./2.*sin(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) + 0.
    rigb1 = -2. * ad1*bd1/tan(a1) + .5*sin(2.*b1)*gd1**2         + 0.
    rigg1 = -2. * ad1*gd1/tan(a1) - 2.*bd1*gd1/tan(b1)           + 0.

    return np.array(
        [ad1,bd1,gd1,
        riga1,rigb1,rigg1])

# Setup system jacobian (Needed since method is implicit)
def dynjac_s3simgeo(state_vec, params):
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,m = params

    #---------- particle 1

    da1riga1 = cos(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + 0.
    db1riga1 = cos(b1)*sin(b1)*sin(2.*a1)*gd1**2.           + 0.
    dg1riga1 = 0.                                           + 0.

    dad1riga1 = 0.                          + 0.
    dbd1riga1 = sin(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sin(2.*a1)*gd1 + 0.

    #----------

    da1rigb1 = 2./sin(a1)**2.*(ad1*bd1) + 0.
    db1rigb1 = cos(2.*b1)*gd1**2.        + 0.
    dg1rigb1 = 0.                        + 0.

    dad1rigb1 = -2./tan(a1)*bd1 + 0.
    dbd1rigb1 = -2./tan(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1    + 0.

    #------------

    da1rigg1 = 2./sin(a1)**2.*(ad1*gd1) + 0.
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)  + 0.
    dg1rigg1 = 0.                        + 0.

    dad1rigg1 = -2./tan(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                  + 0.
    dgd1rigg1 = -2.*(ad1/tan(a1) + bd1/tan(b1)) + 0.


    return np.array([
        [0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, dad1riga1,dbd1riga1,dgd1riga1],
        [da1rigb1,db1rigb1,dg1rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1]
    ])



#######################
# S3 Exact Spring Bar #
#######################

def dynfunc_s3exactbar(state_vec,params):
    l,dl = state_vec
    v,ks,x,m = params
    return np.array([dl,((sin(x/2.)**4. * v**2.)/(sin(l)**2. * tan(l)) - (2.*ks*(l - x/2.))/m)])

# Setup system jacobian (Needed since method is implicit)
def dynjac_s3exactbar(state_vec, params):
    l,dl = state_vec
    v,ks,x,m = params
    return np.array([
        [0.,1.],
        [-((2.*ks)/m) - v**2. * (2. + cos(2.*l)) * 1./sin(l)**4. * sin(x/2.)**4., 0.]
    ])

#####################
# S3 Sim Spring Bar #
#####################

def dynfunc_s3simbar(state_vec,params):
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
    v,ks,x,m = params

    # Distance Function
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # First derivatives of distance function
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)

    spa1 = -(ks*(arccos(d12) - x)*da1d12)/(m*sqrt(1. - d12**2.))
    spb1 = -(ks*(arccos(d12) - x)*db1d12)/(m*sin(a1)**2. * sqrt(1. - d12**2.))
    spg1 = -(ks*(arccos(d12) - x)*dg1d12)/(m*sin(a1)**2. * sin(b1)**2. * sqrt(1. - d12**2.))
    spa2 = -(ks*(arccos(d12) - x)*da2d12)/(m*sqrt(1. - d12**2.))
    spb2 = -(ks*(arccos(d12) - x)*db2d12)/(m*sin(a2)**2. * sqrt(1. - d12**2.))
    spg2 = -(ks*(arccos(d12) - x)*dg2d12)/(m*sin(a2)**2. * sin(b2)**2. * sqrt(1. - d12**2.))

    riga1 = 1./2.*sin(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) - spa1
    rigb1 = -2. * ad1*bd1/tan(a1) + .5*sin(2.*b1)*gd1**2 - spb1
    rigg1 = -2. * ad1*gd1/tan(a1) - 2.*bd1*gd1/tan(b1) - spg1
    riga2 = 1./2.*sin(2. * a2)*(bd2**2. + gd2**2. * sin(b2)**2.) - spa2
    rigb2 = -2. * ad2*bd2/tan(a2) + .5*sin(2.*b2)*gd2**2 - spb2
    rigg2 = -2. * ad2*gd2/tan(a2) - 2.*bd2*gd2/tan(b2) - spg2

    return np.array(
        [ad1,bd1,gd1,ad2,bd2,gd2,
        riga1,rigb1,rigg1,riga2,rigb2,rigg2])

# Setup system jacobian (Needed since method is implicit)
def dynjac_s3simbar(state_vec, params):
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

    def da1D12a1(a1, b1, g1, a2, b2, g2):
        return -cos(a1)*cos(a2) - sin(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

    def db1D12a1(a1, b1, g1, a2, b2, g2):
        return cos(a1)*sin(a2)*(-cos(b2)*sin(b1) + cos(b1)*cos(g1 - g2)*sin(b2))

    def dg1D12a1(a1, b1, g1, a2, b2, g2):
        return -cos(a1)*sin(a2)*sin(b1)*sin(b2)*sin(g1 - g2)

    def da2D12a1(a1, b1, g1, a2, b2, g2):
        return sin(a1)*sin(a2) + cos(a1)*cos(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

    def db2D12a1(a1, b1, g1, a2, b2, g2):
        return cos(a1)*sin(a2)*(cos(b2)*cos(g1 - g2)*sin(b1) - cos(b1)*sin(b2))

    def dg2D12a1(a1, b1, g1, a2, b2, g2):
        return cos(a1)*sin(a2)*sin(b1)*sin(b2)*sin(g1 - g2)

    def db1D12b1(a1, b1, g1, a2, b2, g2):
        return -sin(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

    def dg1D12b1(a1, b1, g1, a2, b2, g2):
        return -cos(b1)*sin(a1)*sin(a2)*sin(b2)*sin(g1 - g2)

    def db2D12b1(a1, b1, g1, a2, b2, g2):
        return sin(a1)*sin(a2)*(cos(b1)*cos(b2)*cos(g1 - g2) + sin(b1)*sin(b2))

    def dg2D12b1(a1, b1, g1, a2, b2, g2):
        return cos(b1)*sin(a1)*sin(a2)*sin(b2)*sin(g1 - g2)

    def dg1D12g1(a1, b1, g1, a2, b2, g2):
        return -cos(g1 - g2)*sin(a1)*sin(a2)*sin(b1)*sin(b2)

    def dg2D12g1(a1, b1, g1, a2, b2, g2):
        return cos(g1 - g2)*sin(a1)*sin(a2)*sin(b1)*sin(b2)
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

    def da2da1V12(m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1):
        # negative sign here so that the term can be added to geo terms later
        return -k/(m*f*sqrt( 1. - d12**2. ))*( (da1d12*da2d12)/sqrt( 1. - d12**2.) + ( arccos(d12) - l )*( da2d12da1 - da1d12*(da2f/f + d12*da2d12/(1. - d12**2.)) ) )
    
    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    v,ks,x,m = params

    # Distance function
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # First derivatives of distance function
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    # Second derivatives of distance function
    da1d12a1=da1D12a1(a1, b1, g1, a2, b2, g2)
    db1d12a1=db1D12a1(a1, b1, g1, a2, b2, g2)
    dg1d12a1=dg1D12a1(a1, b1, g1, a2, b2, g2)
    da2d12a1=da2D12a1(a1, b1, g1, a2, b2, g2)
    db2d12a1=db2D12a1(a1, b1, g1, a2, b2, g2)
    dg2d12a1=dg2D12a1(a1, b1, g1, a2, b2, g2)
    
    da1d12b1=db1d12a1
    db1d12b1=db1D12b1(a1, b1, g1, a2, b2, g2)
    dg1d12b1=dg1D12b1(a1, b1, g1, a2, b2, g2)
    da2d12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12b1=db2D12b1(a1, b1, g1, a2, b2, g2)
    dg2d12b1=dg2D12b1(a1, b1, g1, a2, b2, g2)

    da1d12g1=dg1d12a1
    db1d12g1=dg1d12b1
    dg1d12g1=dg1D12g1(a1, b1, g1, a2, b2, g2)
    da2d12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12g1=dg2D12g1(a1, b1, g1, a2, b2, g2)

    da1d12a2=da2d12a1
    db1d12a2=da2d12b1
    dg1d12a2=da2d12g1
    da2d12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
    db2d12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
    dg2d12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

    da1d12b2=db2d12a1
    db1d12b2=db2d12b1
    dg1d12b2=db2d12g1
    da2d12b2=db2d12a2
    db2d12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

    da1d12g2=dg2d12a1
    db1d12g2=dg2d12b1
    dg1d12g2=dg2d12g1
    da2d12g2=dg2d12a2
    db2d12g2=dg2d12b2
    dg2d12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)

    #---------- particle 1

    da1riga1 = cos(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + da2da1V12(m, 1., ks, x, d12, da1d12, da1d12, 0., da1d12a1)
    db1riga1 = cos(b1)*sin(b1)*sin(2.*a1)*gd1**2.           + da2da1V12(m, 1., ks, x, d12, da1d12, db1d12, 0., db1d12a1)
    dg1riga1 = 0.                                           + da2da1V12(m, 1., ks, x, d12, da1d12, dg1d12, 0., dg1d12a1)

    da2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, da2d12, 0., da2d12a1)
    db2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, db2d12, 0., db2d12a1)
    dg2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, dg2d12, 0., dg2d12a1)

    dad1riga1 = 0.                         + 0.
    dbd1riga1 = sin(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sin(2.*a1)*gd1 + 0.

    dad2riga1 = 0. + 0.
    dbd2riga1 = 0. + 0.
    dgd2riga1 = 0. + 0.

    #----------

    da1rigb1 = 2./sin(a1)**2.*(ad1*bd1) + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, da1d12, sin(2.*a1), da1d12b1)
    db1rigb1 = cos(2.*b1)*gd1**2.       + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, db1d12, 0.,         db1d12b1)
    dg1rigb1 = 0.                       + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, dg1d12, 0.,         dg1d12b1)

    da2rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, da2d12, 0., da2d12b1)
    db2rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, db2d12, 0., db2d12b1)
    dg2rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, dg2d12, 0., dg2d12b1)

    dad1rigb1 = -2./tan(a1)*bd1 + 0.
    dbd1rigb1 = -2./tan(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1   + 0.

    dad2rigb1 = 0. + 0.
    dbd2rigb1 = 0. + 0.
    dgd2rigb1 = 0. + 0.

    #------------

    da1rigg1 = 2./sin(a1)**2.*(ad1*gd1) + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da1d12, sin(2.*a1)*sin(b1)**2., da1d12g1)
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1) + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db1d12, sin(2.*b1)*sin(a1)**2., db1d12g1)
    dg1rigg1 = 0.                       + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg1d12, 0.,                     dg1d12g1)

    da2rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da2d12, 0., da2d12g1)
    db2rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db2d12, 0., db2d12g1)
    dg2rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg2d12, 0., dg2d12g1)

    dad1rigg1 = -2./tan(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                 + 0.
    dgd1rigg1 = -2.*(ad1/tan(a1) + bd1/tan(b1)) + 0.

    dad2rigg1 = 0. + 0.
    dbd2rigg1 = 0. + 0.
    dgd2rigg1 = 0. + 0.

    #------------- particle 2

    da1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, da1d12, 0., da1d12a2)
    db1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, db1d12, 0., db1d12a2)
    dg1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, dg1d12, 0., dg1d12a2)

    da2riga2 = cos(2.*a2) * (bd2**2. + sin(b2)**2.*gd2**2.) + da2da1V12(m, 1., ks, x, d12, da2d12, da2d12, 0., da2d12a2)
    db2riga2 = cos(b2)*sin(b2)*sin(2.*a2)*gd2**2.           + da2da1V12(m, 1., ks, x, d12, da2d12, db2d12, 0., db2d12a2)
    dg2riga2 = 0.                                           + da2da1V12(m, 1., ks, x, d12, da2d12, dg2d12, 0., dg2d12a2)

    dad1riga2 = 0. + 0.
    dbd1riga2 = 0. + 0.
    dgd1riga2 = 0. + 0.

    dad2riga2 = 0.                         + 0.
    dbd2riga2 = sin(2.*a2)*bd2             + 0.
    dgd2riga2 = sin(b2)**2.*sin(2.*a2)*gd2 + 0.

    #----------

    da1rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, da1d12, 0., da1d12b2)
    db1rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, db1d12, 0., db1d12b2)
    dg1rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, dg1d12, 0., dg1d12b2)

    da2rigb2 = 2./sin(a2)**2.*(ad2*bd2) + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, da2d12, sin(2.*a2), da2d12b2)
    db2rigb2 = cos(2.*b2)*gd2**2.       + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, db2d12, 0.,         db2d12b2)
    dg2rigb2 = 0.                       + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, dg2d12, 0.,         dg2d12b2)

    dad1rigb2 = -2./tan(a2)*bd2 + 0.
    dbd1rigb2 = -2./tan(a2)*ad2 + 0.
    dgd1rigb2 = sin(2.*b2)*g2   + 0.

    dad2rigb2 = 0. + 0.
    dbd2rigb2 = 0. + 0.
    dgd2rigb2 = 0. + 0.

    #------------

    da1rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da1d12, 0., da1d12g2)
    db1rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db1d12, 0., db1d12g2)
    dg1rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg1d12, 0., dg1d12g2)

    da2rigg2 = 2./sin(a2)**2.*(ad2*gd2) + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da2d12, sin(2.*a2)*sin(b2)**2., da2d12g2)
    db2rigg2 = 2./sin(b2)**2.*(bd2*gd2) + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db2d12, sin(2.*b2)*sin(a2)**2., db2d12g2)
    dg2rigg2 = 0.                       + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg2d12, 0.,                     dg2d12g2)

    dad1rigg2 = -2./tan(a2)*gd2                 + 0.
    dbd1rigg2 = -2./tan(b2)*gd2                 + 0.
    dgd1rigg2 = -2.*(ad2/tan(a2) + bd2/tan(b2)) + 0.

    dad2rigg2 = 0. + 0.
    dbd2rigg2 = 0. + 0.
    dgd2rigg2 = 0. + 0.

    return np.array([
        [0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, da2riga1,db2riga1,dg2riga1, dad1riga1,dbd1riga1,dgd1riga1, dad2riga1,dbd2riga1,dgd2riga1],
        [da1rigb1,db1rigb1,dg1rigb1, da2rigb1,db2rigb1,dg2rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1, dad2rigb1,dbd2rigb1,dgd2rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, da2rigg1,db2rigg1,dg2rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1, dad2rigg1,dbd2rigg1,dgd2rigg1],
        [da1riga2,db1riga2,dg1riga2, da2riga2,db2riga2,dg2riga2, dad1riga2,dbd1riga2,dgd1riga2, dad2riga2,dbd2riga2,dgd2riga2],
        [da1rigb2,db1rigb2,dg1rigb2, da2rigb2,db2rigb2,dg2rigb2, dad1rigb2,dbd1rigb2,dgd1rigb2, dad2rigb2,dbd2rigb2,dgd2rigb2],
        [da1rigg2,db1rigg2,dg1rigg2, da2rigg2,db2rigg2,dg2rigg2, dad1rigg2,dbd1rigg2,dgd1rigg2, dad2rigg2,dbd2rigg2,dgd2rigg2]
    ])


############################
# S3 Exact Spring Triangle #
############################

def dynfunc_s3exacttriangle(state_vec,params):
    l,H,dl,dH = state_vec
    v,ks,x,m = params
    P = m*v*(1. + 2.*cos(x/2.)**2.)
    aux1 = cos(l)*cos(H)
    aux2 = sin(l)*sin(H)
    return np.array([
        dl,dH,
        -(((2.*ks*(-x + arccos(aux1))*cos(H)*sin(l))/sqrt(1. - aux1**2.) + (2.*(-dH*m + P)**2.*cos(l)*sin(l))/(m*(1. + 2.*cos(l)**2.)**2.) + (4.*ks*(-x + 2.*l)*cos(l)*sin(l))/sqrt(1. - cos(2.*l)**2.))/(2.*m)),
        ((2.*ks*(-x + arccos(aux1))*cos(l)*sin(H))/sqrt(1. - aux1**2.) + (2.*dl*(-dH*m + P)*sin(2.*l))/(2. + cos(2.*l))**2.)/(m*(-1. + 1./(2. + cos(2.*l))))])

# Setup system jacobian (Needed since method is implicit)
def dynjac_s3exacttriangle(state_vec, params):
    l,H,dl,dH = state_vec
    v,ks,x,m = params
    P = m*v*(1. + 2.*cos(x/2.)**2.)
    aux1 = cos(l)*cos(H)
    aux2 = sin(l)*sin(H)
    return np.array([
        [0.,0.,1.,0.],
        [0.,0.,0.,1.],
        [-(1./(2.*m))*((2.*(-dH*m + P)**2.*cos(2.*l))/(m*(1. + 2.*cos(l)**2.)**2.) + (2.*ks*(-x + arccos(aux1))*aux1)/sqrt(1. - aux1**2.) + (16.*(-dH*m + P)**2.*cos(l)**2.*sin(l)**2.)/(m*(1. + 2.*cos(l)**2.)**3.) - (2.*ks*(-x + arccos(aux1))*cos(H)**3.*cos(l)*sin(l)**2.)/(1. - aux1**2. )**(3./2.) + (2.*ks*cos(H)**2.*sin(l)**2.)/(1. - aux1**2.) - (16.*ks*(-x + 2.*l)*cos(l)**2.*sin(l)**2.*(cos(2.*l)))/(1. - (cos(2.*l))**2.)**(3./2.) + (16.*ks*cos(l)**2.*sin(l)**2.)/(1. - (cos(2.*l))**2.) + (4.*ks*(-x + 2.*l)*cos(2.*l))/sqrt(1. - (cos(2.*l))**2.)),
        -((ks*(x - arccos(aux1) + aux1*sqrt(1. - aux1**2.))*aux2)/(m*(1. - aux1**2.)**(3./2.))),
         0.,
         ((-dH*m + P)*sin(2.*l))/(m*(2. + cos(2.*l))**2.)],
        [-((2.*m*sin(2.*l)*((2.*ks*(-x + arccos(aux1))*cos(l)*sin(H))/sqrt(1. - aux1**2.) + (2.*dl*(-dH*m + P)*sin(2.*l))/(2. + cos(2.*l))**2.))/((2. + cos(2.*l))**2.*(-m + m/(2. + cos(2.*l)))**2.)) + (1./(-m + m/(2. + cos(2.*l))))*((4.*dl*(-dH*m + P)*cos(2.*l))/(2. + cos(2.*l))**2. - (2.*ks*(-x + arccos(aux1))*aux1**2.*aux2)/(1. - aux1**2.)**(3./2.) + (2.*ks*aux1*aux2)/(1. - aux1**2.) - (2.*ks*(-x + arccos(aux1))*aux2)/sqrt(1. - aux1**2.) + (8.*dl*(-dH*m + P)*sin(2.*l)**2.)/(2. + cos(2.*l))**3.),
         -((ks*(2. + cos(2.*l))*(sqrt(1. - aux1**2.)*sin(H)**2. + (-x + arccos(aux1))*cos(H)*sin(l)*tan(l)))/(m*(1. - aux1**2.)**(3./2.))),
         (2.*(dH*m - P)*tan(l))/(m*(2. + cos(2.*l))),
         (2.*dl*tan(l))/(2. + cos(2.*l))]
        ])



##########################
# S3 Sim Spring Triangle #
##########################

def dynfunc_s3simtriangle(state_vec,params):
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
    
    a1,b1,g1,a2,b2,g2,a3,b3,g3,ad1,bd1,gd1,ad2,bd2,gd2,ad3,bd3,gd3 = state_vec
    v,ks,x,m = params

    # Distance Functions
    # spring 12
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # spring 13
    d13 = D12(a1, b1, g1, a3, b3, g3)
    # spring 23
    d23 = D12(a2, b2, g2, a3, b3, g3)
    # First derivatives of distance function
    # spring 12
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    # spring 13
    da1d13 = da1D12(a1, b1, g1, a3, b3, g3)
    db1d13 = db1D12(a1, b1, g1, a3, b3, g3)
    dg1d13 = dg1D12(a1, b1, g1, a3, b3, g3)
    da3d13 = da1D12(a3, b3, g3, a1, b1, g1)
    db3d13 = db1D12(a3, b3, g3, a1, b1, g1)
    dg3d13 = dg1D12(a3, b3, g3, a1, b1, g1)
    # spring 23
    da2d23 = da1D12(a2, b2, g2, a3, b3, g3)
    db2d23 = db1D12(a2, b2, g2, a3, b3, g3)
    dg2d23 = dg1D12(a2, b2, g2, a3, b3, g3)
    da3d23 = da1D12(a3, b3, g3, a2, b2, g2)
    db3d23 = db1D12(a3, b3, g3, a2, b2, g2)
    dg3d23 = dg1D12(a3, b3, g3, a2, b2, g2)

    def da1V12(m, f, k, l, d12, da1d12):
        # negative sign here so that the term can be added to geo terms later (but positive due to double negative)
        return (k*(arccos(d12) - l)*da1d12)/(m*f*sqrt( 1. - d12**2. ))

    
    spa1 = da1V12(m,1.,ks,x,d12,da1d12)                        + da1V12(m,1.,ks,x,d13,da1d13)
    spb1 = da1V12(m,sin(a1)**2.,ks,x,d12,db1d12)               + da1V12(m,sin(a1)**2.,ks,x,d13,db1d13)
    spg1 = da1V12(m,sin(a1)**2.*sin(b1)**2.,ks,x,d12,dg1d12)   + da1V12(m,sin(a1)**2.*sin(b1)**2.,ks,x,d13,dg1d13)
    spa2 = da1V12(m,1.,ks,x,d12,da2d12)                        + da1V12(m,1.,ks,x,d23,da2d23)
    spb2 = da1V12(m,sin(a2)**2.,ks,x,d12,db2d12)               + da1V12(m,sin(a2)**2.,ks,x,d23,db2d23)
    spg2 = da1V12(m,sin(a2)**2.*sin(b2)**2.,ks,x,d12,dg2d12)   + da1V12(m,sin(a2)**2.*sin(b2)**2.,ks,x,d23,dg2d23)
    spa3 = da1V12(m,1.,ks,x,d13,da3d13)                        + da1V12(m,1.,ks,x,d23,da3d23)
    spb3 = da1V12(m,sin(a3)**2.,ks,x,d13,db3d13)               + da1V12(m,sin(a3)**2.,ks,x,d23,db3d23)
    spg3 = da1V12(m,sin(a3)**2.*sin(b3)**2.,ks,x,d13,dg3d13)   + da1V12(m,sin(a3)**2.*sin(b3)**2.,ks,x,d23,dg3d23)

    riga1 = 1./2.*sin(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.)   + spa1
    rigb1 = -2. * ad1*bd1/tan(a1) + .5*sin(2.*b1)*gd1**2           + spb1
    rigg1 = -2. * ad1*gd1/tan(a1) - 2.*bd1*gd1/tan(b1)             + spg1
    riga2 = 1./2.*sin(2. * a2)*(bd2**2. + gd2**2. * sin(b2)**2.)   + spa2
    rigb2 = -2. * ad2*bd2/tan(a2) + .5*sin(2.*b2)*gd2**2           + spb2
    rigg2 = -2. * ad2*gd2/tan(a2) - 2.*bd2*gd2/tan(b2)             + spg2
    riga3 = 1./2.*sin(2. * a3)*(bd3**2. + gd3**2. * sin(b3)**2.)   + spa3
    rigb3 = -2. * ad3*bd3/tan(a3) + .5*sin(2.*b3)*gd3**2           + spb3
    rigg3 = -2. * ad3*gd3/tan(a3) - 2.*bd3*gd3/tan(b3)             + spg3

    return np.array(
        [ad1,bd1,gd1,ad2,bd2,gd2,ad3,bd3,gd3,
        riga1,rigb1,rigg1,riga2,rigb2,rigg2,riga3,rigb3,rigg3])

# Setup system jacobian (Needed since method is implicit)
def dynjac_s3simtriangle(state_vec, params):
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

    def da1D12a1(a1, b1, g1, a2, b2, g2):
        return -cos(a1)*cos(a2) - sin(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

    def db1D12a1(a1, b1, g1, a2, b2, g2):
        return cos(a1)*sin(a2)*(-cos(b2)*sin(b1) + cos(b1)*cos(g1 - g2)*sin(b2))

    def dg1D12a1(a1, b1, g1, a2, b2, g2):
        return -cos(a1)*sin(a2)*sin(b1)*sin(b2)*sin(g1 - g2)

    def da2D12a1(a1, b1, g1, a2, b2, g2):
        return sin(a1)*sin(a2) + cos(a1)*cos(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

    def db2D12a1(a1, b1, g1, a2, b2, g2):
        return cos(a1)*sin(a2)*(cos(b2)*cos(g1 - g2)*sin(b1) - cos(b1)*sin(b2))

    def dg2D12a1(a1, b1, g1, a2, b2, g2):
        return cos(a1)*sin(a2)*sin(b1)*sin(b2)*sin(g1 - g2)

    def db1D12b1(a1, b1, g1, a2, b2, g2):
        return -sin(a1)*sin(a2)*(cos(b1)*cos(b2) + cos(g1 - g2)*sin(b1)*sin(b2))

    def dg1D12b1(a1, b1, g1, a2, b2, g2):
        return -cos(b1)*sin(a1)*sin(a2)*sin(b2)*sin(g1 - g2)

    def db2D12b1(a1, b1, g1, a2, b2, g2):
        return sin(a1)*sin(a2)*(cos(b1)*cos(b2)*cos(g1 - g2) + sin(b1)*sin(b2))

    def dg2D12b1(a1, b1, g1, a2, b2, g2):
        return cos(b1)*sin(a1)*sin(a2)*sin(b2)*sin(g1 - g2)

    def dg1D12g1(a1, b1, g1, a2, b2, g2):
        return -cos(g1 - g2)*sin(a1)*sin(a2)*sin(b1)*sin(b2)

    def dg2D12g1(a1, b1, g1, a2, b2, g2):
        return cos(g1 - g2)*sin(a1)*sin(a2)*sin(b1)*sin(b2)
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

    def da2da1V12(m, f, k, l, d12, da1d12, da2d12, da2f, da2d12da1):
        # negative sign here so that the term can be added to geo terms later (positive due to double negative)
        return k/(m*f*sqrt( 1. - d12**2. ))*( (da1d12*da2d12)/sqrt( 1. - d12**2.) + ( arccos(d12) - l )*( da2d12da1 - da1d12*(da2f/f + d12*da2d12/(1. - d12**2.)) ) )

    a1,b1,g1,a2,b2,g2,a3,b3,g3,ad1,bd1,gd1,ad2,bd2,gd2,ad3,bd3,gd3 = state_vec
    v,ks,x,m = params

    # Distance Functions
    # spring 12
    d12 = D12(a1, b1, g1, a2, b2, g2)
    # spring 13
    d13 = D12(a1, b1, g1, a3, b3, g3)
    # spring 23
    d23 = D12(a2, b2, g2, a3, b3, g3)
    # First derivatives of distance function
    # spring 12
    da1d12 = da1D12(a1, b1, g1, a2, b2, g2)
    db1d12 = db1D12(a1, b1, g1, a2, b2, g2)
    dg1d12 = dg1D12(a1, b1, g1, a2, b2, g2)
    da2d12 = da1D12(a2, b2, g2, a1, b1, g1)
    db2d12 = db1D12(a2, b2, g2, a1, b1, g1)
    dg2d12 = dg1D12(a2, b2, g2, a1, b1, g1)
    # spring 13
    da1d13 = da1D12(a1, b1, g1, a3, b3, g3)
    db1d13 = db1D12(a1, b1, g1, a3, b3, g3)
    dg1d13 = dg1D12(a1, b1, g1, a3, b3, g3)
    da3d13 = da1D12(a3, b3, g3, a1, b1, g1)
    db3d13 = db1D12(a3, b3, g3, a1, b1, g1)
    dg3d13 = dg1D12(a3, b3, g3, a1, b1, g1)
    # spring 23
    da2d23 = da1D12(a2, b2, g2, a3, b3, g3)
    db2d23 = db1D12(a2, b2, g2, a3, b3, g3)
    dg2d23 = dg1D12(a2, b2, g2, a3, b3, g3)
    da3d23 = da1D12(a3, b3, g3, a2, b2, g2)
    db3d23 = db1D12(a3, b3, g3, a2, b2, g2)
    dg3d23 = dg1D12(a3, b3, g3, a2, b2, g2)
    # Second derivatives of distance function
    # spring 12
    da1d12a1=da1D12a1(a1, b1, g1, a2, b2, g2)
    db1d12a1=db1D12a1(a1, b1, g1, a2, b2, g2)
    dg1d12a1=dg1D12a1(a1, b1, g1, a2, b2, g2)
    da2d12a1=da2D12a1(a1, b1, g1, a2, b2, g2)
    db2d12a1=db2D12a1(a1, b1, g1, a2, b2, g2)
    dg2d12a1=dg2D12a1(a1, b1, g1, a2, b2, g2)
    
    da1d12b1=db1d12a1
    db1d12b1=db1D12b1(a1, b1, g1, a2, b2, g2)
    dg1d12b1=dg1D12b1(a1, b1, g1, a2, b2, g2)
    da2d12b1 = db2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12b1=db2D12b1(a1, b1, g1, a2, b2, g2)
    dg2d12b1=dg2D12b1(a1, b1, g1, a2, b2, g2)

    da1d12g1=dg1d12a1
    db1d12g1=dg1d12b1
    dg1d12g1=dg1D12g1(a1, b1, g1, a2, b2, g2)
    da2d12g1 = dg2D12a1(a2, b2, g2, a1, b1, g1)
    db2d12g1 = dg2D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12g1=dg2D12g1(a1, b1, g1, a2, b2, g2)

    da1d12a2=da2d12a1
    db1d12a2=da2d12b1
    dg1d12a2=da2d12g1
    da2d12a2 = da1D12a1(a2, b2, g2, a1, b1, g1)
    db2d12a2 = db1D12a1(a2, b2, g2, a1, b1, g1)
    dg2d12a2 = dg1D12a1(a2, b2, g2, a1, b1, g1)

    da1d12b2=db2d12a1
    db1d12b2=db2d12b1
    dg1d12b2=db2d12g1
    da2d12b2=db2d12a2
    db2d12b2 = db1D12b1(a2, b2, g2, a1, b1, g1)
    dg2d12b2 = dg1D12b1(a2, b2, g2, a1, b1, g1)

    da1d12g2=dg2d12a1
    db1d12g2=dg2d12b1
    dg1d12g2=dg2d12g1
    da2d12g2=dg2d12a2
    db2d12g2=dg2d12b2
    dg2d12g2 = dg1D12g1(a2, b2, g2, a1, b1, g1)
    # spring 13
    da1d13a1=da1D12a1(a1, b1, g1, a3, b3, g3)
    db1d13a1=db1D12a1(a1, b1, g1, a3, b3, g3)
    dg1d13a1=dg1D12a1(a1, b1, g1, a3, b3, g3)
    da3d13a1=da2D12a1(a1, b1, g1, a3, b3, g3)
    db3d13a1=db2D12a1(a1, b1, g1, a3, b3, g3)
    dg3d13a1=dg2D12a1(a1, b1, g1, a3, b3, g3)
    
    da1d13b1=db1d13a1
    db1d13b1=db1D12b1(a1, b1, g1, a3, b3, g3)
    dg1d13b1=dg1D12b1(a1, b1, g1, a3, b3, g3)
    da3d13b1 = db2D12a1(a3, b3, g3, a1, b1, g1)
    db3d13b1=db2D12b1(a1, b1, g1, a3, b3, g3)
    dg3d13b1=dg2D12b1(a1, b1, g1, a3, b3, g3)

    da1d13g1=dg1d13a1
    db1d13g1=dg1d13b1
    dg1d13g1=dg1D12g1(a1, b1, g1, a3, b3, g3)
    da3d13g1 = dg2D12a1(a3, b3, g3, a1, b1, g1)
    db3d13g1 = dg2D12b1(a3, b3, g3, a1, b1, g1)
    dg3d13g1=dg2D12g1(a1, b1, g1, a3, b3, g3)

    da1d13a3=da3d13a1
    db1d13a3=da3d13b1
    dg1d13a3=da3d13g1
    da3d13a3 = da1D12a1(a3, b3, g3, a1, b1, g1)
    db3d13a3 = db1D12a1(a3, b3, g3, a1, b1, g1)
    dg3d13a3 = dg1D12a1(a3, b3, g3, a1, b1, g1)

    da1d13b3=db3d13a1
    db1d13b3=db3d13b1
    dg1d13b3=db3d13g1
    da3d13b3=db3d13a3
    db3d13b3 = db1D12b1(a3, b3, g3, a1, b1, g1)
    dg3d13b3 = dg1D12b1(a3, b3, g3, a1, b1, g1)

    da1d13g3=dg3d13a1
    db1d13g3=dg3d13b1
    dg1d13g3=dg3d13g1
    da3d13g3=dg3d13a3
    db3d13g3=dg3d13b3
    dg3d13g3 = dg1D12g1(a3, b3, g3, a1, b1, g1)
    # spring 23
    da2d23a2=da1D12a1(a2, b2, g2, a3, b3, g3)
    db2d23a2=db1D12a1(a2, b2, g2, a3, b3, g3)
    dg2d23a2=dg1D12a1(a2, b2, g2, a3, b3, g3)
    da3d23a2=da2D12a1(a2, b2, g2, a3, b3, g3)
    db3d23a2=db2D12a1(a2, b2, g2, a3, b3, g3)
    dg3d23a2=dg2D12a1(a2, b2, g2, a3, b3, g3)

    da2d23b2=db2d23a2
    db2d23b2=db1D12b1(a2, b2, g2, a3, b3, g3)
    dg2d23b2=dg1D12b1(a2, b2, g2, a3, b3, g3)
    da3d23b2 = db2D12a1(a3, b3, g3, a2, b2, g2)
    db3d23b2=db2D12b1(a2, b2, g2, a3, b3, g3)
    dg3d23b2=dg2D12b1(a2, b2, g2, a3, b3, g3)

    da2d23g2=dg2d23a2
    db2d23g2=dg2d23b2
    dg2d23g2=dg1D12g1(a2, b2, g2, a3, b3, g3)
    da3d23g2 = dg2D12a1(a3, b3, g3, a2, b2, g2)
    db3d23g2 = dg2D12b1(a3, b3, g3, a2, b2, g2)
    dg3d23g2=dg2D12g1(a2, b2, g2, a3, b3, g3)

    da2d23a3=da3d23a2
    db2d23a3=da3d23b2
    dg2d23a3=da3d23g2
    da3d23a3 = da1D12a1(a3, b3, g3, a2, b2, g2)
    db3d23a3 = db1D12a1(a3, b3, g3, a2, b2, g2)
    dg3d23a3 = dg1D12a1(a3, b3, g3, a2, b2, g2)

    da2d23b3=db3d23a2
    db2d23b3=db3d23b2
    dg2d23b3=db3d23g2
    da3d23b3=db3d23a3
    db3d23b3 = db1D12b1(a3, b3, g3, a2, b2, g2)
    dg3d23b3 = dg1D12b1(a3, b3, g3, a2, b2, g2)

    da2d23g3=dg3d23a2
    db2d23g3=dg3d23b2
    dg2d23g3=dg3d23g2
    da3d23g3=dg3d23a3
    db3d23g3=dg3d23b3
    dg3d23g3 = dg1D12g1(a3, b3, g3, a2, b2, g2)

    #---------- particle 1

    da1riga1 = cos(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + da2da1V12(m, 1., ks, x, d12, da1d12, da1d12, 0., da1d12a1) + da2da1V12(m, 1., ks, x, d13, da1d13, da1d13, 0., da1d13a1)
    db1riga1 = cos(b1)*sin(b1)*sin(2.*a1)*gd1**2.           + da2da1V12(m, 1., ks, x, d12, da1d12, db1d12, 0., db1d12a1) + da2da1V12(m, 1., ks, x, d13, da1d13, db1d13, 0., db1d13a1)
    dg1riga1 = 0.                                           + da2da1V12(m, 1., ks, x, d12, da1d12, dg1d12, 0., dg1d12a1) + da2da1V12(m, 1., ks, x, d13, da1d13, dg1d13, 0., dg1d13a1)

    da2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, da2d12, 0., da2d12a1)
    db2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, db2d12, 0., db2d12a1)
    dg2riga1 = 0. + da2da1V12(m, 1., ks, x, d12, da1d12, dg2d12, 0., dg2d12a1)

    da3riga1 = 0. + da2da1V12(m, 1., ks, x, d13, da1d12, da3d13, 0., da3d13a1)
    db3riga1 = 0. + da2da1V12(m, 1., ks, x, d13, da1d12, db3d13, 0., db3d13a1)
    dg3riga1 = 0. + da2da1V12(m, 1., ks, x, d13, da1d12, dg3d13, 0., dg3d13a1)

    dad1riga1 = 0.                         + 0.
    dbd1riga1 = sin(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sin(2.*a1)*gd1 + 0.

    dad2riga1 = 0. + 0.
    dbd2riga1 = 0. + 0.
    dgd2riga1 = 0. + 0.

    dad3riga1 = 0. + 0.
    dbd3riga1 = 0. + 0.
    dgd3riga1 = 0. + 0.

    #----------

    da1rigb1 = 2./sin(a1)**2.*(ad1*bd1) + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, da1d12, sin(2.*a1), da1d12b1) + da2da1V12(m, sin(a1)**2., ks, x, d13, db1d13, da1d13, sin(2.*a1), da1d13b1)
    db1rigb1 = cos(2.*b1)*gd1**2.       + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, db1d12, 0.,         db1d12b1) + da2da1V12(m, sin(a1)**2., ks, x, d13, db1d13, db1d13, 0.,         db1d13b1)
    dg1rigb1 = 0.                       + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, dg1d12, 0.,         dg1d12b1) + da2da1V12(m, sin(a1)**2., ks, x, d13, db1d13, dg1d13, 0.,         dg1d13b1)

    da2rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, da2d12, 0., da2d12b1)
    db2rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, db2d12, 0., db2d12b1)
    dg2rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d12, db1d12, dg2d12, 0., dg2d12b1)

    da3rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d13, db1d13, da3d13, 0., da3d13b1)
    db3rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d13, db1d13, db3d13, 0., db3d13b1)
    dg3rigb1 = 0. + da2da1V12(m, sin(a1)**2., ks, x, d13, db1d13, dg3d13, 0., dg3d13b1)

    dad1rigb1 = -2./tan(a1)*bd1 + 0.
    dbd1rigb1 = -2./tan(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1    + 0.

    dad2rigb1 = 0. + 0.
    dbd2rigb1 = 0. + 0.
    dgd2rigb1 = 0. + 0.

    dad3rigb1 = 0. + 0.
    dbd3rigb1 = 0. + 0.
    dgd3rigb1 = 0. + 0.

    #------------

    da1rigg1 = 2./sin(a1)**2.*(ad1*gd1) + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da1d12, sin(2.*a1)*sin(b1)**2., da1d12g1) + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, da1d13, sin(2.*a1)*sin(b1)**2., da1d13g1)
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1) + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db1d12, sin(2.*b1)*sin(a1)**2., db1d12g1) + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, db1d13, sin(2.*b1)*sin(a1)**2., db1d13g1)
    dg1rigg1 = 0.                       + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg1d12, 0.,                     dg1d12g1) + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, dg1d13, 0.,                     dg1d13g1)

    da2rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, da2d12, 0., da2d12g1)
    db2rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, db2d12, 0., db2d12g1)
    dg2rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d12, dg1d12, dg2d12, 0., dg2d12g1)

    da3rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, da3d13, 0., da3d13g1)
    db3rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, db3d13, 0., db3d13g1)
    dg3rigg1 = 0. + da2da1V12(m, sin(a1)**2.*sin(b1)**2., ks, x, d13, dg1d13, dg3d13, 0., dg3d13g1)

    dad1rigg1 = -2./tan(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                 + 0.
    dgd1rigg1 = -2.*(ad1/tan(a1) + bd1/tan(b1)) + 0.

    dad2rigg1 = 0. + 0.
    dbd2rigg1 = 0. + 0.
    dgd2rigg1 = 0. + 0.

    dad3rigg1 = 0. + 0.
    dbd3rigg1 = 0. + 0.
    dgd3rigg1 = 0. + 0.

    #------------- particle 2

    da1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, da1d12, 0., da1d12a2)
    db1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, db1d12, 0., db1d12a2)
    dg1riga2 = 0. + da2da1V12(m, 1., ks, x, d12, da2d12, dg1d12, 0., dg1d12a2)

    da2riga2 = cos(2.*a2) * (bd2**2. + sin(b2)**2.*gd2**2.) + da2da1V12(m, 1., ks, x, d12, da2d12, da2d12, 0., da2d12a2) + da2da1V12(m, 1., ks, x, d23, da2d23, da2d23, 0., da2d23a2)
    db2riga2 = cos(b2)*sin(b2)*sin(2.*a2)*gd2**2.           + da2da1V12(m, 1., ks, x, d12, da2d12, db2d12, 0., db2d12a2) + da2da1V12(m, 1., ks, x, d23, da2d23, db2d23, 0., db2d23a2)
    dg2riga2 = 0.                                           + da2da1V12(m, 1., ks, x, d12, da2d12, dg2d12, 0., dg2d12a2) + da2da1V12(m, 1., ks, x, d23, da2d23, dg2d23, 0., dg2d23a2)

    da3riga2 = 0. + da2da1V12(m, 1., ks, x, d23, da2d23, da3d23, 0., da3d23a2)
    db3riga2 = 0. + da2da1V12(m, 1., ks, x, d23, da2d23, db3d23, 0., db3d23a2)
    dg3riga2 = 0. + da2da1V12(m, 1., ks, x, d23, da2d23, dg3d23, 0., dg3d23a2)
    
    dad1riga2 = 0. + 0.
    dbd1riga2 = 0. + 0.
    dgd1riga2 = 0. + 0.

    dad2riga2 = 0.                         + 0.
    dbd2riga2 = sin(2.*a2)*bd2             + 0.
    dgd2riga2 = sin(b2)**2.*sin(2.*a2)*gd2 + 0.

    dad3riga2 = 0. + 0.
    dbd3riga2 = 0. + 0.
    dgd3riga2 = 0. + 0.

    #----------

    da1rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, da1d12, 0., da1d12b2)
    db1rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, db1d12, 0., db1d12b2)
    dg1rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, dg1d12, 0., dg1d12b2)

    da2rigb2 = 2./sin(a2)**2.*(ad2*bd2) + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, da2d12, sin(2.*a2), da2d12b2) + da2da1V12(m, sin(a2)**2., ks, x, d23, db2d23, da2d23, sin(2.*a2), da2d23b2)
    db2rigb2 = cos(2.*b2)*gd2**2.       + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, db2d12, 0.,         db2d12b2) + da2da1V12(m, sin(a2)**2., ks, x, d23, db2d23, db2d23, 0.,         db2d23b2)
    dg2rigb2 = 0.                       + da2da1V12(m, sin(a2)**2., ks, x, d12, db2d12, dg2d12, 0.,         dg2d12b2) + da2da1V12(m, sin(a2)**2., ks, x, d23, db2d23, dg2d23, 0.,         dg2d23b2)

    da3rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d23, db2d23, da3d23, 0., da3d23b2)
    db3rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d23, db2d23, db3d23, 0., db3d23b2)
    dg3rigb2 = 0. + da2da1V12(m, sin(a2)**2., ks, x, d23, db2d23, dg3d23, 0., dg3d23b2)
    
    dad1rigb2 = -2./tan(a2)*bd2 + 0.
    dbd1rigb2 = -2./tan(a2)*ad2 + 0.
    dgd1rigb2 = sin(2.*b2)*g2   + 0.

    dad2rigb2 = 0. + 0.
    dbd2rigb2 = 0. + 0.
    dgd2rigb2 = 0. + 0.

    dad3rigb2 = 0. + 0.
    dbd3rigb2 = 0. + 0.
    dgd3rigb2 = 0. + 0.

    #------------

    da1rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da1d12, 0., da1d12g2)
    db1rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db1d12, 0., db1d12g2)
    dg1rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg1d12, 0., dg1d12g2)

    da2rigg2 = 2./sin(a2)**2.*(ad2*gd2) + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, da2d12, sin(2.*a2)*sin(b2)**2., da2d12g2) + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, da2d23, sin(2.*a2)*sin(b2)**2., da2d23g2)
    db2rigg2 = 2./sin(b2)**2.*(bd2*gd2) + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, db2d12, sin(2.*b2)*sin(a2)**2., db2d12g2) + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, db2d23, sin(2.*b2)*sin(a2)**2., db2d23g2)
    dg2rigg2 = 0.                       + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d12, dg2d12, dg2d12, 0.,                     dg2d12g2) + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, dg2d23, 0.,                     dg2d23g2)

    da3rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, da3d23, 0., da3d23g2)
    db3rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, db3d23, 0., db3d23g2)
    dg3rigg2 = 0. + da2da1V12(m, sin(a2)**2.*sin(b2)**2., ks, x, d23, dg2d23, dg3d23, 0., dg3d23g2)
    
    dad1rigg2 = -2./tan(a2)*gd2                 + 0.
    dbd1rigg2 = -2./tan(b2)*gd2                 + 0.
    dgd1rigg2 = -2.*(ad2/tan(a2) + bd2/tan(b2)) + 0.

    dad2rigg2 = 0. + 0.
    dbd2rigg2 = 0. + 0.
    dgd2rigg2 = 0. + 0.

    dad3rigg2 = 0. + 0.
    dbd3rigg2 = 0. + 0.
    dgd3rigg2 = 0. + 0.

    #------------- particle 3

    da1riga3 = 0. + da2da1V12(m, 1., ks, x, d13, da3d13, da1d13, 0., da1d13a3)
    db1riga3 = 0. + da2da1V12(m, 1., ks, x, d13, da3d13, db1d13, 0., db1d13a3)
    dg1riga3 = 0. + da2da1V12(m, 1., ks, x, d13, da3d13, dg1d13, 0., dg1d13a3)

    da2riga3 = 0. + da2da1V12(m, 1., ks, x, d23, da3d23, da2d23, 0., da2d23a3)
    db2riga3 = 0. + da2da1V12(m, 1., ks, x, d23, da3d23, db2d23, 0., db2d23a3)
    dg2riga3 = 0. + da2da1V12(m, 1., ks, x, d23, da3d23, dg2d23, 0., dg2d23a3)

    da3riga3 = cos(2.*a3) * (bd3**2. + sin(b3)**2.*gd3**2.) + da2da1V12(m, 1., ks, x, d12, da3d13, da3d13, 0., da3d13a3) + da2da1V12(m, 1., ks, x, d23, da3d23, da3d23, 0., da3d23a3)
    db3riga3 = cos(b3)*sin(b3)*sin(2.*a3)*gd3**2.           + da2da1V12(m, 1., ks, x, d12, da3d13, db3d13, 0., db3d13a3) + da2da1V12(m, 1., ks, x, d23, da3d23, db3d23, 0., db3d23a3)
    dg3riga3 = 0.                                           + da2da1V12(m, 1., ks, x, d12, da3d13, dg3d13, 0., dg3d13a3) + da2da1V12(m, 1., ks, x, d23, da3d23, dg3d23, 0., dg3d23a3)
    
    dad1riga3 = 0. + 0.
    dbd1riga3 = 0. + 0.
    dgd1riga3 = 0. + 0.

    dad2riga3 = 0. + 0.
    dbd2riga3 = 0. + 0.
    dgd2riga3 = 0. + 0.

    dad3riga3 = 0.                         + 0.
    dbd3riga3 = sin(2.*a3)*bd3             + 0.
    dgd3riga3 = sin(b3)**2.*sin(2.*a3)*gd3 + 0.

    #----------

    da1rigb3 = 0. + da2da1V12(m, sin(a3)**2., ks, x, d13, db3d13, da1d13, 0., da1d13b3)
    db1rigb3 = 0. + da2da1V12(m, sin(a3)**2., ks, x, d13, db3d13, db1d13, 0., db1d13b3)
    dg1rigb3 = 0. + da2da1V12(m, sin(a3)**2., ks, x, d13, db3d13, dg1d13, 0., dg1d13b3)

    da2rigb3 = 0. + da2da1V12(m, sin(a3)**2., ks, x, d23, db3d23, da2d23, 0., da2d23b3)
    db2rigb3 = 0. + da2da1V12(m, sin(a3)**2., ks, x, d23, db3d23, db2d23, 0., db2d23b3)
    dg2rigb3 = 0. + da2da1V12(m, sin(a3)**2., ks, x, d23, db3d23, dg2d23, 0., dg2d23b3)

    da3rigb3 = 2./sin(a3)**2.*(ad3*bd3) + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, da3d13, sin(2.*a3), da3d13b3) + da2da1V12(m, sin(a3)**2., ks, x, d23, db3d23, da3d23, sin(2.*a3), da3d23b3)
    db3rigb3 = cos(2.*b3)*gd3**2.       + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, db3d13, 0.,         db3d13b3) + da2da1V12(m, sin(a3)**2., ks, x, d23, db3d23, db3d23, 0.,         db3d23b3)
    dg3rigb3 = 0.                       + da2da1V12(m, sinh(a3)**2., ks, x, d13, db3d13, dg3d13, 0.,         dg3d13b3) + da2da1V12(m, sin(a3)**2., ks, x, d23, db3d23, dg3d23, 0.,         dg3d23b3)

    dad1rigb3 = 0. + 0.
    dbd1rigb3 = 0. + 0.
    dgd1rigb3 = 0. + 0.

    dad2rigb3 = 0. + 0.
    dbd2rigb3 = 0. + 0.
    dgd2rigb3 = 0. + 0.

    dad3rigb3 = -2./tan(a3)*bd3 + 0.
    dbd3rigb3 = -2./tan(a3)*ad3 + 0.
    dgd3rigb3 = sin(2.*b3)*g3    + 0.

    #------------

    da1rigg3 = 0. + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, da1d13, 0., da1d13g3)
    db1rigg3 = 0. + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, db1d13, 0., db1d13g3)
    dg1rigg3 = 0. + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, dg1d13, 0., dg1d13g3)

    da2rigg3 = 0. + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, da2d23, 0., da2d23g3)
    db2rigg3 = 0. + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, db2d23, 0., db2d23g3)
    dg2rigg3 = 0. + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, dg2d23, 0., dg2d23g3)

    da3rigg3 = 2./sin(a3)**2.*(ad3*gd3) + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, da3d13, sin(2.*a3)*sin(b3)**2., da3d13g3) + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, da3d23, sin(2.*a3)*sin(b3)**2., da3d23g3)
    db3rigg3 = 2./sin(b3)**2.*(bd3*gd3) + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, db3d13, sin(2.*b3)*sin(a3)**2., db3d13g3) + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, db3d23, sin(2.*b3)*sin(a3)**2., db3d23g3)
    dg3rigg3 = 0.                       + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d13, dg3d13, dg3d13, 0.,                     dg3d13g3) + da2da1V12(m, sin(a3)**2.*sin(b3)**2., ks, x, d23, dg3d23, dg3d23, 0.,                     dg3d23g3)

    dad1rigg3 = 0. + 0.
    dbd1rigg3 = 0. + 0.
    dgd1rigg3 = 0. + 0.

    dad2rigg3 = 0. + 0.
    dbd2rigg3 = 0. + 0.
    dgd2rigg3 = 0. + 0.

    dad3rigg3 = -2./tan(a3)*gd3                 + 0.
    dbd3rigg3 = -2./tan(b3)*gd3                 + 0.
    dgd3rigg3 = -2.*(ad3/tan(a3) + bd3/tan(b3)) + 0.

    return np.array([
        [0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.],
        [da1riga1,db1riga1,dg1riga1, da2riga1,db2riga1,dg2riga1, da3riga1,db3riga1,dg3riga1, dad1riga1,dbd1riga1,dgd1riga1, dad2riga1,dbd2riga1,dgd2riga1, dad3riga1,dbd3riga1,dgd3riga1],
        [da1rigb1,db1rigb1,dg1rigb1, da2rigb1,db2rigb1,dg2rigb1, da3rigb1,db3rigb1,dg3rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1, dad2rigb1,dbd2rigb1,dgd2rigb1, dad3rigb1,dbd3rigb1,dgd3rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, da2rigg1,db2rigg1,dg2rigg1, da3rigg1,db3rigg1,dg3rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1, dad2rigg1,dbd2rigg1,dgd2rigg1, dad3rigg1,dbd3rigg1,dgd3rigg1],
        [da1riga2,db1riga2,dg1riga2, da2riga2,db2riga2,dg2riga2, da3riga2,db3riga2,dg3riga2, dad1riga2,dbd1riga2,dgd1riga2, dad2riga2,dbd2riga2,dgd2riga2, dad3riga2,dbd3riga2,dgd3riga2],
        [da1rigb2,db1rigb2,dg1rigb2, da2rigb2,db2rigb2,dg2rigb2, da3rigb2,db3rigb2,dg3rigb2, dad1rigb2,dbd1rigb2,dgd1rigb2, dad2rigb2,dbd2rigb2,dgd2rigb2, dad3rigb2,dbd3rigb2,dgd3rigb2],
        [da1rigg2,db1rigg2,dg1rigg2, da2rigg2,db2rigg2,dg2rigg2, da3rigg2,db3rigg2,dg3rigg2, dad1rigg2,dbd1rigg2,dgd1rigg2, dad2rigg2,dbd2rigg2,dgd2rigg2, dad3rigg2,dbd3rigg2,dgd3rigg2],
        [da1riga3,db1riga3,dg1riga3, da2riga3,db2riga3,dg2riga3, da3riga3,db3riga3,dg3riga3, dad1riga3,dbd1riga3,dgd1riga3, dad2riga3,dbd2riga3,dgd2riga3, dad3riga3,dbd3riga3,dgd3riga3],
        [da1rigb3,db1rigb3,dg1rigb3, da2rigb3,db2rigb3,dg2rigb3, da3rigb3,db3rigb3,dg3rigb3, dad1rigb3,dbd1rigb3,dgd1rigb3, dad2rigb3,dbd2rigb3,dgd2rigb3, dad3rigb3,dbd3rigb3,dgd3rigb3],
        [da1rigg3,db1rigg3,dg1rigg3, da2rigg3,db2rigg3,dg2rigg3, da3rigg3,db3rigg3,dg3rigg3, dad1rigg3,dbd1rigg3,dgd1rigg3, dad2rigg3,dbd2rigg3,dgd2rigg3, dad3rigg3,dbd3rigg3,dgd3rigg3]
    ])


############################
# S3 Sim Central Potential #
############################

def dynfunc_s3simcenpot(state_vec,params):
    
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,G,ms,m = params

    # Central Potential first derivative
    # Based on potential converging at negative source at antipode
    def s3da1cenpot(state_vec,params):
        a1,b1,g1,ad1,bd1,gd1 = state_vec
        v,G,ms,m = params
        return -G*ms/sin(a1)**2.

    riga1 = 1./2.*sin(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) + s3da1cenpot(state_vec,params)
    rigb1 = -2. * ad1*bd1/tan(a1) + .5*sin(2.*b1)*gd1**2         + 0.
    rigg1 = -2. * ad1*gd1/tan(a1) - 2.*bd1*gd1/tan(b1)           + 0.

    return np.array(
        [ad1,bd1,gd1,
        riga1,rigb1,rigg1])

# Setup system jacobian (Needed since method is implicit)
def dynjac_s3simcenpot(state_vec, params):
    a1,b1,g1,ad1,bd1,gd1 = state_vec
    v,G,ms,m = params

    # Central Potential Second derivative
    def s3da1cenpotda1(state_vec,params):
        a1,b1,g1,ad1,bd1,gd1 = state_vec
        v,G,ms,m = params
        return 2.*G*ms/sin(a1)**2./tan(a1)**2.

    #---------- particle 1

    da1riga1 = cos(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.) + s3da1cenpotda1(state_vec,params)
    db1riga1 = cos(b1)*sin(b1)*sin(2.*a1)*gd1**2.           + 0.
    dg1riga1 = 0.                                           + 0.

    dad1riga1 = 0.                          + 0.
    dbd1riga1 = sin(2.*a1)*bd1             + 0.
    dgd1riga1 = sin(b1)**2.*sin(2.*a1)*gd1 + 0.

    #----------

    da1rigb1 = 2./sin(a1)**2.*(ad1*bd1) + 0.
    db1rigb1 = cos(2.*b1)*gd1**2.        + 0.
    dg1rigb1 = 0.                        + 0.

    dad1rigb1 = -2./tan(a1)*bd1 + 0.
    dbd1rigb1 = -2./tan(a1)*ad1 + 0.
    dgd1rigb1 = sin(2.*b1)*g1    + 0.

    #------------

    da1rigg1 = 2./sin(a1)**2.*(ad1*gd1) + 0.
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)  + 0.
    dg1rigg1 = 0.                        + 0.

    dad1rigg1 = -2./tan(a1)*gd1                 + 0.
    dbd1rigg1 = -2./tan(b1)*gd1                  + 0.
    dgd1rigg1 = -2.*(ad1/tan(a1) + bd1/tan(b1)) + 0.


    return np.array([
        [0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, dad1riga1,dbd1riga1,dgd1riga1],
        [da1rigb1,db1rigb1,dg1rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1]
    ])


#############################
# S3 Sim 2 Sphere Collision #
#############################

def dynfunc_s3sim2ballcol(state_vec,params):
    
    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    m1,m2,r1,r2 = params

    riga1 = 1./2.*sin(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.)
    rigb1 = -2. * ad1*bd1/tan(a1) + .5*sin(2.*b1)*gd1**2
    rigg1 = -2. * ad1*gd1/tan(a1) - 2.*bd1*gd1/tan(b1)

    riga2 = 1./2.*sin(2. * a2)*(bd2**2. + gd2**2. * sin(b2)**2.)
    rigb2 = -2. * ad2*bd2/tan(a2) + .5*sin(2.*b2)*gd2**2
    rigg2 = -2. * ad2*gd2/tan(a2) - 2.*bd2*gd2/tan(b2)

    return np.array(
        [ad1,bd1,gd1,ad2,bd2,gd2,
        riga1,rigb1,rigg1,riga2,rigb2,rigg2])

# Setup system jacobian (Needed since method is implicit)
def dynjac_s3sim2ballcol(state_vec, params):
    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    m1,m2,r1,r2 = params

    da1riga1 = cos(2.*a1) * (bd1**2. + sin(b1)**2.*gd1**2.)
    db1riga1 = cos(b1)*sin(b1)*sin(2.*a1)*gd1**2.          
    dg1riga1 = 0.                                          

    da2riga1 = 0.
    db2riga1 = 0.
    dg2riga1 = 0.

    dad1riga1 = 0.                        
    dbd1riga1 = sin(2.*a1)*bd1            
    dgd1riga1 = sin(b1)**2.*sin(2.*a1)*gd1

    dad2riga1 = 0.
    dbd2riga1 = 0.
    dgd2riga1 = 0.

    #----------

    da1rigb1 = 2./sin(a1)**2.*(ad1*bd1)
    db1rigb1 = cos(2.*b1)*gd1**2.      
    dg1rigb1 = 0.                      

    da2rigb1 = 0.
    db2rigb1 = 0.
    dg2rigb1 = 0.

    dad1rigb1 = -2./tan(a1)*bd1
    dbd1rigb1 = -2./tan(a1)*ad1
    dgd1rigb1 = sin(2.*b1)*g1  

    dad2rigb1 = 0.
    dbd2rigb1 = 0.
    dgd2rigb1 = 0.

    #------------

    da1rigg1 = 2./sin(a1)**2.*(ad1*gd1)
    db1rigg1 = 2./sin(b1)**2.*(bd1*gd1)
    dg1rigg1 = 0.                      

    da2rigg1 = 0.
    db2rigg1 = 0.
    dg2rigg1 = 0.

    dad1rigg1 = -2./tan(a1)*gd1                
    dbd1rigg1 = -2./tan(b1)*gd1                
    dgd1rigg1 = -2.*(ad1/tan(a1) + bd1/tan(b1))

    dad2rigg1 = 0.
    dbd2rigg1 = 0.
    dgd2rigg1 = 0.

    #------------- particle 2

    da1riga2 = 0.
    db1riga2 = 0.
    dg1riga2 = 0.

    da2riga2 = cos(2.*a2) * (bd2**2. + sin(b2)**2.*gd2**2.)
    db2riga2 = cos(b2)*sin(b2)*sin(2.*a2)*gd2**2.          
    dg2riga2 = 0.                                          

    dad1riga2 = 0.
    dbd1riga2 = 0.
    dgd1riga2 = 0.

    dad2riga2 = 0.                        
    dbd2riga2 = sin(2.*a2)*bd2            
    dgd2riga2 = sin(b2)**2.*sin(2.*a2)*gd2

    #----------

    da1rigb2 = 0.
    db1rigb2 = 0.
    dg1rigb2 = 0.

    da2rigb2 = 2./sin(a2)**2.*(ad2*bd2)
    db2rigb2 = cos(2.*b2)*gd2**2.      
    dg2rigb2 = 0.                      

    dad1rigb2 = -2./tan(a2)*bd2
    dbd1rigb2 = -2./tan(a2)*ad2
    dgd1rigb2 = sin(2.*b2)*g2  

    dad2rigb2 = 0.
    dbd2rigb2 = 0.
    dgd2rigb2 = 0.

    #------------

    da1rigg2 = 0.
    db1rigg2 = 0.
    dg1rigg2 = 0.

    da2rigg2 = 2./sin(a2)**2.*(ad2*gd2)
    db2rigg2 = 2./sin(b2)**2.*(bd2*gd2)
    dg2rigg2 = 0.                      

    dad1rigg2 = -2./tan(a2)*gd2                
    dbd1rigg2 = -2./tan(b2)*gd2                
    dgd1rigg2 = -2.*(ad2/tan(a2) + bd2/tan(b2))

    dad2rigg2 = 0.
    dbd2rigg2 = 0.
    dgd2rigg2 = 0.

    return np.array([
        [0., 0., 0.,  0., 0., 0.,  1., 0., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 1., 0.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 1.,  0., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  1., 0., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 1., 0.],
        [0., 0., 0.,  0., 0., 0.,  0., 0., 0.,  0., 0., 1.],
        [da1riga1,db1riga1,dg1riga1, da2riga1,db2riga1,dg2riga1, dad1riga1,dbd1riga1,dgd1riga1, dad2riga1,dbd2riga1,dgd2riga1],
        [da1rigb1,db1rigb1,dg1rigb1, da2rigb1,db2rigb1,dg2rigb1, dad1rigb1,dbd1rigb1,dgd1rigb1, dad2rigb1,dbd2rigb1,dgd2rigb1],
        [da1rigg1,db1rigg1,dg1rigg1, da2rigg1,db2rigg1,dg2rigg1, dad1rigg1,dbd1rigg1,dgd1rigg1, dad2rigg1,dbd2rigg1,dgd2rigg1],
        [da1riga2,db1riga2,dg1riga2, da2riga2,db2riga2,dg2riga2, dad1riga2,dbd1riga2,dgd1riga2, dad2riga2,dbd2riga2,dgd2riga2],
        [da1rigb2,db1rigb2,dg1rigb2, da2rigb2,db2rigb2,dg2rigb2, dad1rigb2,dbd1rigb2,dgd1rigb2, dad2rigb2,dbd2rigb2,dgd2rigb2],
        [da1rigg2,db1rigg2,dg1rigg2, da2rigg2,db2rigg2,dg2rigg2, dad1rigg2,dbd1rigg2,dgd1rigg2, dad2rigg2,dbd2rigg2,dgd2rigg2]
    ])




#################################################################
#################### Bump Space Test Systems ####################
#################################################################


#####################################
# Bump Space Sim Curvature Detector #
#####################################

def dynfunc_bssimcurdet(state_vec,params):
    
    # Conformal Factors from Metric
    def cofac(a1, b1, g1):
        return a1**2. + b1**2. + g1**2. + 2.
    # Conformal Factors Time Derivative
    def dtcofac(a1, b1, g1, ad1, bd1, gd1):
        return 2.*a1*ad1 + 2.*b1*bd1 + 2.*g1*gd1
    

    # PPQ Terms
    def ppqfunc(a1, b1, g1, a2, b2, g2):
        return sqrt((a1 + a2)**2. + (b1 + b2)**2. + (g1 + g2)**2.)
    
    def dtppqfunc(a1, b1, g1, a2, b2, g2, ad1, bd1, gd1, ad2, bd2, gd2):
        return (2.*(a1 + a2)*(ad1 + ad2) + 2.*(b1 + b2)*(bd1 + bd2) + 2.*(g1 + g2)*(gd1 + gd2))/(2.*sqrt((a1 + a2)**2. + (b1 + b2)**2. + (g1 + g2)**2.))
    
    def da1ppqfunc(a1, b1, g1, a2, b2, g2):
        return (2.*(a1 + a2)*(1.))/(2.*sqrt((a1 + a2)**2. + (b1 + b2)**2. + (g1 + g2)**2.))
    
    def db1ppqfunc(a1, b1, g1, a2, b2, g2):
        return (2.*(b1 + b2)*(1.))/(2.*sqrt((a1 + a2)**2. + (b1 + b2)**2. + (g1 + g2)**2.))
    
    def dg1ppqfunc(a1, b1, g1, a2, b2, g2):
        return (2.*(g1 + g2)*(1.))/(2.*sqrt((a1 + a2)**2. + (b1 + b2)**2. + (g1 + g2)**2.))
    
    # For the remaining three functions use:
    # da2ppq = da1ppqfunc(a2, b2, g2, a1, b1, g1)
    # db2ppq = db1ppqfunc(a2, b2, g2, a1, b1, g1)
    # dg2ppq = dg1ppqfunc(a2, b2, g2, a1, b1, g1)
    

    # PNQ Terms
    def pnqfunc(a1, b1, g1, a2, b2, g2):
        return sqrt((a1 - a2)**2. + (b1 - b2)**2. + (g1 - g2)**2.)
    
    def dtpnqfunc(a1, b1, g1, a2, b2, g2):
        return (2.*(a1 - a2)*(ad1 - ad2) + 2.*(b1 - b2)*(bd1 - bd2) + 2.*(g1 - g2)*(gd1 - gd2))/(2.*sqrt((a1 - a2)**2. + (b1 - b2)**2. + (g1 - g2)**2.))
    
    def da1pnqfunc(a1, b1, g1, a2, b2, g2):
        return (2.*(a1 - a2)*(1.))/(2.*sqrt((a1 - a2)**2. + (b1 - b2)**2. + (g1 - g2)**2.))
    
    def db1pnqfunc(a1, b1, g1, a2, b2, g2):
        return (2.*(a1 - a2)*(1.))/(2.*sqrt((a1 - a2)**2. + (b1 - b2)**2. + (g1 - g2)**2.))
    
    def dg1pnqfunc(a1, b1, g1, a2, b2, g2):
        return (2.*(a1 - a2)*(1.))/(2.*sqrt((a1 - a2)**2. + (b1 - b2)**2. + (g1 - g2)**2.))
    
    # For the remaining three functions use:
    # da2pnq = da1pnqfunc(a2, b2, g2, a1, b1, g1)
    # db2pnq = db1pnqfunc(a2, b2, g2, a1, b1, g1)
    # dg2pnq = dg1pnqfunc(a2, b2, g2, a1, b1, g1)

    
    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec
    ks,x,m = params

    # PPQ & PNQ Functions
    ppq = ppqfunc(a1, b1, g1, a2, b2, g2)
    pnq = pnqfunc(a1, b1, g1, a2, b2, g2)
    # First derivatives of PPQ & PNQ functions
    da1ppq = da1ppqfunc(a1, b1, g1, a2, b2, g2)
    db1ppq = db1ppqfunc(a1, b1, g1, a2, b2, g2)
    dg1ppq = dg1ppqfunc(a1, b1, g1, a2, b2, g2)
    da2ppq = da1ppqfunc(a2, b2, g2, a1, b1, g1)
    db2ppq = db1ppqfunc(a2, b2, g2, a1, b1, g1)
    dg2ppq = dg1ppqfunc(a2, b2, g2, a1, b1, g1)

    da1pnq = da1pnqfunc(a1, b1, g1, a2, b2, g2)
    db1pnq = db1pnqfunc(a1, b1, g1, a2, b2, g2)
    dg1pnq = dg1pnqfunc(a1, b1, g1, a2, b2, g2)
    da2pnq = da1pnqfunc(a2, b2, g2, a1, b1, g1)
    db2pnq = db1pnqfunc(a2, b2, g2, a1, b1, g1)
    dg2pnq = dg1pnqfunc(a2, b2, g2, a1, b1, g1)


    spa1 = -(ks*(arccos(d12) - x)*da1d12)/(m*sqrt(1. - d12**2.))
    spb1 = -(ks*(arccos(d12) - x)*db1d12)/(m*sin(a1)**2. * sqrt(1. - d12**2.))
    spg1 = -(ks*(arccos(d12) - x)*dg1d12)/(m*sin(a1)**2. * sin(b1)**2. * sqrt(1. - d12**2.))
    spa2 = -(ks*(arccos(d12) - x)*da2d12)/(m*sqrt(1. - d12**2.))
    spb2 = -(ks*(arccos(d12) - x)*db2d12)/(m*sin(a2)**2. * sqrt(1. - d12**2.))
    spg2 = -(ks*(arccos(d12) - x)*dg2d12)/(m*sin(a2)**2. * sin(b2)**2. * sqrt(1. - d12**2.))

    riga1 = 1./2.*sin(2. * a1)*(bd1**2. + gd1**2. * sin(b1)**2.) - spa1
    rigb1 = -2. * ad1*bd1/tan(a1) + .5*sin(2.*b1)*gd1**2 - spb1
    rigg1 = -2. * ad1*gd1/tan(a1) - 2.*bd1*gd1/tan(b1) - spg1
    riga2 = 1./2.*sin(2. * a2)*(bd2**2. + gd2**2. * sin(b2)**2.) - spa2
    rigb2 = -2. * ad2*bd2/tan(a2) + .5*sin(2.*b2)*gd2**2 - spb2
    rigg2 = -2. * ad2*gd2/tan(a2) - 2.*bd2*gd2/tan(b2) - spg2

    return np.array(
        [ad1,bd1,gd1,ad2,bd2,gd2,
        riga1,rigb1,rigg1,riga2,rigb2,rigg2])









