from sympy import *
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3exactbar, dynjac_h3exactbar, dynfunc_h3simbar, dynjac_h3simbar, dynfunc_s3exactbar, dynjac_s3exactbar, dynfunc_s3simbar, dynjac_s3simbar


a1n,b1n,g1n,a2n,b2n,g2n = sp.symbols('a1n b1n g1n a2n b2n g2n')
ad1n,bd1n,gd1n,ad2n,bd2n,gd2n = sp.symbols('ad1n bd1n gd1n ad2n bd2n gd2n')
a1c1,a1c2,b1c1,b1c2,g1c1,g1c2,a2c1,a2c2,b2c1,b2c2,g2c1,g2c2 = sp.symbols('a1c1 a1c2 b1c1 b1c2 g1c1 g1c2 a2c1 a2c2 b2c1 b2c2 g2c1 g2c2')
ad1c1,ad1c2,bd1c1,bd1c2,gd1c1,gd1c2,ad2c1,ad2c2,bd2c1,bd2c2,gd2c1,gd2c2 = sp.symbols('ad1c1 ad1c2 bd1c1 bd1c2 gd1c1 gd1c2 ad2c1 ad2c2 bd2c1 bd2c2 gd2c1 gd2c2')
m1,m2,leq,dt,a11,a12,a21,a22,lamc1,lamc2 = sp.symbols('m1 m2 leq dt a11 a12 a21 a22 lamc1 lamc2')

d12c1 = cosh(a1c1)*cosh(a2c1) - sinh(a1c1)*cos(b1c1)*sinh(a2c1)*cos(b2c1) - sinh(a1c1)*sin(b1c1)*sinh(a2c1)*sin(b2c1)*cos(g1c1 - g2c1)
d12c2 = cosh(a1c2)*cosh(a2c2) - sinh(a1c2)*cos(b1c2)*sinh(a2c2)*cos(b2c2) - sinh(a1c2)*sin(b1c2)*sinh(a2c2)*sin(b2c2)*cos(g1c2 - g2c2)

dtd12c2 = cosh(a2c2)*(ad1c2 - ad2c2*(cos(b1c2)*cos(b2c2) + cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2)))*sinh(a1c2) + (cosh(a1c2)*(ad2c2 - ad1c2*(cos(b1c2)*cos(b2c2) + cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2))) + (cos(b2c2)*(bd1c2 - bd2c2*cos(g1c2 - g2c2))*sin(b1c2) + sin(b2c2)*(cos(b1c2)*(bd2c2 - bd1c2*cos(g1c2 - g2c2)) + (gd1c2 - gd2c2)*sin(b1c2)*sin(g1c2 - g2c2)))*sinh(a1c2))*sinh(a2c2)

a1termc1 = 1./2.*sinh(2. * a1c1)*(bd1c1**2. + gd1c1**2. * sin(b1c1)**2.)
a1termc2 = 1./2.*sinh(2. * a1c2)*(bd1c2**2. + gd1c2**2. * sin(b1c2)**2.)
a2termc1 = 1./2.*sinh(2. * a2c1)*(bd2c1**2. + gd2c1**2. * sin(b2c1)**2.)
a2termc2 = 1./2.*sinh(2. * a2c2)*(bd2c2**2. + gd2c2**2. * sin(b2c2)**2.)
b1termc1 = -2. * ad1c1*bd1c1/tanh(a1c1) + .5*sin(2.*b1c1)*gd1c1**2
b1termc2 = -2. * ad1c2*bd1c2/tanh(a1c2) + .5*sin(2.*b1c2)*gd1c2**2
b2termc1 = -2. * ad2c1*bd2c1/tanh(a2c1) + .5*sin(2.*b2c1)*gd2c1**2
b2termc2 = -2. * ad2c2*bd2c2/tanh(a2c2) + .5*sin(2.*b2c2)*gd2c2**2
g1termc1 = -2. * ad1c1*gd1c1/tanh(a1c1) - 2.*bd1c1*gd1c1/tan(b1c1)
g1termc2 = -2. * ad1c2*gd1c2/tanh(a1c2) - 2.*bd1c2*gd1c2/tan(b1c2)
g2termc1 = -2. * ad2c1*gd2c1/tanh(a2c1) - 2.*bd2c1*gd2c1/tan(b2c1)
g2termc2 = -2. * ad2c2*gd2c2/tanh(a2c2) - 2.*bd2c2*gd2c2/tan(b2c2)

#dtdtd12c2 = ad1c2**2.*cosh(a1c2)*cosh(a2c2) + ad2c2**2.*cosh(a1c2)*cosh(a2c2) - 2.*ad1c2*ad2c2*cosh(a1c2)*cosh(a2c2)*(cos(b1c2)*cos(b2c2) + cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2)) - 2.*ad2c2*cosh(a2c2)*(-bd1c2*cos(b2c2)*sin(b1c2) + bd2c2*cos(b2c2)*cos(g1c2 - g2c2)*sin(b1c2) - bd2c2*cos(b1c2)*sin(b2c2) + bd1c2*cos(b1c2)*cos(g1c2 - g2c2)*sin(b2c2) - (gd1c2 - gd2c2)*sin(b1c2)*sin(b2c2)*sin(g1c2 - g2c2))*sinh(a1c2) - 2.*ad1c2*cosh(a1c2)*(-bd1c2*cos(b2c2)*sin(b1c2) + bd2c2*cos(b2c2)*cos(g1c2 - g2c2)*sin(b1c2) - bd2c2*cos(b1c2)*sin(b2c2) + bd1c2*cos(b1c2)*cos(g1c2 - g2c2)*sin(b2c2) - (gd1c2 - gd2c2)*sin(b1c2)*sin(b2c2)*sin(g1c2 - g2c2))*sinh(a2c2) + 2.*ad1c2*ad2c2*sinh(a1c2)*sinh(a2c2) - ad1c2**2.*(cos(b1c2)*cos(b2c2) + cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2))*sinh(a1c2)*sinh(a2c2) - ad2c2**2.*(cos(b1c2)*cos(b2c2) + cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2))*sinh(a1c2)*sinh(a2c2) + cosh(a2c2)*sinh[a1c2]*(a1termc2 - cona1c2) - cosh(a1c2)*(cos(b1c2)*cos(b2c2) + cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2))*sinh(a2c2)*(a1termc2 - cona1c2) - cosh(a2c2)*(cos(b1c2)*cos(b2c2) + cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2))*sinh(a1c2)*(a2termc2 - cona2c2) + cosh(a1c2)*sinh(a2c2)*(a2termc2 - cona2c2) - sinh(a1c2)*sinh(a2c2)*(-bd1c2**2.*cos(b1c2)*cos(b2c2) - bd2c2**2.*cos(b1c2)*cos(b2c2) + 2.*bd1c2*bd2c2*cos(b1c2)*cos(b2c2)*cos(g1c2 - g2c2) + 2.*bd1c2*bd2c2*sin(b1c2)*sin(b2c2) - bd1c2**2.*cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2) - bd2c2**2.*cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2) - (gd1c2 - gd2c2)**2.*cos(g1c2 - g2c2)*sin(b1c2)*sin(b2c2) - 2.*bd2c2*(gd1c2 - gd2c2)*cos(b2c2)*sin(b1c2)*sin(g1c2 - g2c2) - 2.*bd1c2*(gd1c2 - gd2c2)*cos(b1c2)*sin(b2c2)*sin(g1c2 - g2c2) - cos(b2c2)*sin(b1c2)*(b1termc2 - conb1c2) + cos(b1c2)*cos(g1c2 - g2c2)*sin(b2c2)*(b1termc2 - conb1c2) + cos(b2c2)*cos(g1c2 - g2c2)*sin(b1c2)*(b2termc2 - conb2c2) - cos(b1c2)*sin(b2c2)*(b2termc2 - conb2c2) - sin(b1c2)*sin(b2c2)*sin(g1c2 - g2c2)*((g1termc2 - cong1c2) - (g2termc2 - cong2c2)))

cona1c1 = (lamc1*d12c1.diff(a1c1))/(m1*sqrt(d12c1**2. - 1))
cona1c2 = (lamc2*d12c2.diff(a1c2))/(m1*sqrt(d12c2**2. - 1))
cona2c1 = (lamc1*d12c1.diff(a2c1))/(m2*sqrt(d12c1**2. - 1))
cona2c2 = (lamc2*d12c2.diff(a2c2))/(m2*sqrt(d12c2**2. - 1))
conb1c1 = (lamc1*d12c1.diff(b1c1))/(m1*sinh(a1c1)**2.*sqrt(d12c1**2. - 1))
conb1c2 = (lamc2*d12c2.diff(b1c2))/(m1*sinh(a1c2)**2.*sqrt(d12c2**2. - 1))
conb2c1 = (lamc1*d12c1.diff(b2c1))/(m2*sinh(a2c1)**2.*sqrt(d12c1**2. - 1))
conb2c2 = (lamc2*d12c2.diff(b2c2))/(m2*sinh(a2c2)**2.*sqrt(d12c2**2. - 1))
cong1c1 = (lamc1*d12c1.diff(g1c1))/(m1*sinh(a1c1)**2.*sin(b1c1)**2.*sqrt(d12c1**2. - 1))
cong1c2 = (lamc2*d12c2.diff(g1c2))/(m1*sinh(a1c2)**2.*sin(b1c2)**2.*sqrt(d12c2**2. - 1))
cong2c1 = (lamc1*d12c1.diff(g2c1))/(m2*sinh(a2c1)**2.*sin(b2c1)**2.*sqrt(d12c1**2. - 1))
cong2c2 = (lamc2*d12c2.diff(g2c2))/(m2*sinh(a2c2)**2.*sin(b2c2)**2.*sqrt(d12c2**2. - 1))

ex1 = a1c1 - a1n - dt*(a11*ad1c1 + a12*ad1c2)
ex2 = a1c2 - a1n - dt*(a21*ad1c1 + a22*ad1c2)

ex3 = b1c1 - b1n - dt*(a11*bd1c1 + a12*bd1c2)
ex4 = b1c2 - b1n - dt*(a21*bd1c1 + a22*bd1c2)

ex5 = g1c1 - g1n - dt*(a11*gd1c1 + a12*gd1c2)
ex6 = g1c2 - g1n - dt*(a21*gd1c1 + a22*gd1c2)

ex7 = a2c1 - a2n - dt*(a11*ad2c1 + a12*ad2c2)
ex8 = a2c2 - a2n - dt*(a21*ad2c1 + a22*ad2c2)

ex9 = b2c1 - b2n - dt*(a11*bd2c1 + a12*bd2c2)
ex10 = b2c2 - b2n - dt*(a21*bd2c1 + a22*bd2c2)

ex11 = g2c1 - g2n - dt*(a11*gd2c1 + a12*gd2c2)
ex12 = g2c2 - g2n - dt*(a21*gd2c1 + a22*gd2c2)


ex13 = ad1c1 - ad1n - dt*(a11*(a1termc1 - cona1c1) + a12*(a1termc2 - cona1c2))
ex14 = ad1c2 - ad1n - dt*(a21*(a1termc1 - cona1c1) + a22*(a1termc2 - cona1c2))

ex15 = bd1c1 - bd1n - dt*(a11*(b1termc1 - conb1c1) + a12*(b1termc2 - conb1c2))
ex16 = bd1c2 - bd1n - dt*(a21*(b1termc1 - conb1c1) + a22*(b1termc2 - conb1c2))

ex17 = gd1c1 - gd1n - dt*(a11*(g1termc1 - cong1c1) + a12*(g1termc2 - cong1c2))
ex18 = gd1c2 - gd1n - dt*(a21*(g1termc1 - cong1c1) + a22*(g1termc2 - cong1c2))

ex19 = ad2c1 - ad2n - dt*(a11*(a2termc1 - cona2c1) + a12*(a2termc2 - cona2c2))
ex20 = ad2c2 - ad2n - dt*(a21*(a2termc1 - cona2c1) + a22*(a2termc2 - cona2c2))

ex21 = bd2c1 - bd2n - dt*(a11*(b2termc1 - conb2c1) + a12*(b2termc2 - conb2c2))
ex22 = bd2c2 - bd2n - dt*(a21*(b2termc1 - conb2c1) + a22*(b2termc2 - conb2c2))

ex23 = gd2c1 - gd2n - dt*(a11*(g2termc1 - cong2c1) + a12*(g2termc2 - cong2c2))
ex24 = gd2c2 - gd2n - dt*(a21*(g2termc1 - cong2c1) + a22*(g2termc2 - cong2c2))

ex25 = acosh(d12c2) - leq
ex26 = dtd12c2/sqrt(d12c2**2. - 1)

inputvals={a1n:0,a1c1:0,a1c2:0,b1n:0,b1c1:0,b1c2:0,g1n:0,g1c1:0,g1c2:0,a2n:0,a2c1:0,a2c2:0,b2n:0,b2c1:0,b2c2:0,g2n:0,g2c1:0,g2c2:0,ad1n:0,ad1c1:0,ad1c2:0,bd1n:0,bd1c1:0,bd1c2:0,gd1n:0,gd1c1:0,gd1c2:0,ad2n:0,ad2c1:0,ad2c2:0,bd2n:0,bd2c1:0,bd2c2:0,gd2n:0,gd2c1:0,gd2c2:0,lamc1:0,lamc2:0,m1:1,m2:1,dt:1,leq:1,a11:5./12.,a12:-1./12.,a21:3./4,a22:1./4.}


def evaljac(arr,inputval):
    boo = arr.copy()
    row,col = boo.shape
    for a in range(row):
        for b in range(col):
            print(a,b)
            boo[a,b] = boo[a,b].evalf(subs=inputval)
    return boo

def evalcon(arr,inputval):
    boo = arr.copy()
    row = len(boo)
    for a in range(row):
            print(a)
            boo[a] = boo[a].evalf(subs=inputval)
    return boo

jacobian = np.array([
        [ex1.diff(a1c1),
         ex1.diff(a1c2),
         ex1.diff(b1c1),
         ex1.diff(b1c2),
         ex1.diff(g1c1),
         ex1.diff(g1c2),
         ex1.diff(a2c1),
         ex1.diff(a2c2),
         ex1.diff(b2c1),
         ex1.diff(b2c2),
         ex1.diff(g2c1),
         ex1.diff(g2c2),
         ex1.diff(ad1c1),
         ex1.diff(ad1c2),
         ex1.diff(bd1c1),
         ex1.diff(bd1c2),
         ex1.diff(gd1c1),
         ex1.diff(gd1c2),
         ex1.diff(ad2c1),
         ex1.diff(ad2c2),
         ex1.diff(bd2c1),
         ex1.diff(bd2c2),
         ex1.diff(gd2c1),
         ex1.diff(gd2c2),
         ex1.diff(lamc1),
         ex1.diff(lamc2)],
        [ex2.diff(a1c1),
         ex2.diff(a1c2),
         ex2.diff(b1c1),
         ex2.diff(b1c2),
         ex2.diff(g1c1),
         ex2.diff(g1c2),
         ex2.diff(a2c1),
         ex2.diff(a2c2),
         ex2.diff(b2c1),
         ex2.diff(b2c2),
         ex2.diff(g2c1),
         ex2.diff(g2c2),
         ex2.diff(ad1c1),
         ex2.diff(ad1c2),
         ex2.diff(bd1c1),
         ex2.diff(bd1c2),
         ex2.diff(gd1c1),
         ex2.diff(gd1c2),
         ex2.diff(ad2c1),
         ex2.diff(ad2c2),
         ex2.diff(bd2c1),
         ex2.diff(bd2c2),
         ex2.diff(gd2c1),
         ex2.diff(gd2c2),
         ex2.diff(lamc1),
         ex2.diff(lamc2)],
        [ex3.diff(a1c1),
         ex3.diff(a1c2),
         ex3.diff(b1c1),
         ex3.diff(b1c2),
         ex3.diff(g1c1),
         ex3.diff(g1c2),
         ex3.diff(a2c1),
         ex3.diff(a2c2),
         ex3.diff(b2c1),
         ex3.diff(b2c2),
         ex3.diff(g2c1),
         ex3.diff(g2c2),
         ex3.diff(ad1c1),
         ex3.diff(ad1c2),
         ex3.diff(bd1c1),
         ex3.diff(bd1c2),
         ex3.diff(gd1c1),
         ex3.diff(gd1c2),
         ex3.diff(ad2c1),
         ex3.diff(ad2c2),
         ex3.diff(bd2c1),
         ex3.diff(bd2c2),
         ex3.diff(gd2c1),
         ex3.diff(gd2c2),
         ex3.diff(lamc1),
         ex3.diff(lamc2)],
        [ex4.diff(a1c1),
         ex4.diff(a1c2),
         ex4.diff(b1c1),
         ex4.diff(b1c2),
         ex4.diff(g1c1),
         ex4.diff(g1c2),
         ex4.diff(a2c1),
         ex4.diff(a2c2),
         ex4.diff(b2c1),
         ex4.diff(b2c2),
         ex4.diff(g2c1),
         ex4.diff(g2c2),
         ex4.diff(ad1c1),
         ex4.diff(ad1c2),
         ex4.diff(bd1c1),
         ex4.diff(bd1c2),
         ex4.diff(gd1c1),
         ex4.diff(gd1c2),
         ex4.diff(ad2c1),
         ex4.diff(ad2c2),
         ex4.diff(bd2c1),
         ex4.diff(bd2c2),
         ex4.diff(gd2c1),
         ex4.diff(gd2c2),
         ex4.diff(lamc1),
         ex4.diff(lamc2)],
        [ex5.diff(a1c1),
         ex5.diff(a1c2),
         ex5.diff(b1c1),
         ex5.diff(b1c2),
         ex5.diff(g1c1),
         ex5.diff(g1c2),
         ex5.diff(a2c1),
         ex5.diff(a2c2),
         ex5.diff(b2c1),
         ex5.diff(b2c2),
         ex5.diff(g2c1),
         ex5.diff(g2c2),
         ex5.diff(ad1c1),
         ex5.diff(ad1c2),
         ex5.diff(bd1c1),
         ex5.diff(bd1c2),
         ex5.diff(gd1c1),
         ex5.diff(gd1c2),
         ex5.diff(ad2c1),
         ex5.diff(ad2c2),
         ex5.diff(bd2c1),
         ex5.diff(bd2c2),
         ex5.diff(gd2c1),
         ex5.diff(gd2c2),
         ex5.diff(lamc1),
         ex5.diff(lamc2)],
        [ex6.diff(a1c1),
         ex6.diff(a1c2),
         ex6.diff(b1c1),
         ex6.diff(b1c2),
         ex6.diff(g1c1),
         ex6.diff(g1c2),
         ex6.diff(a2c1),
         ex6.diff(a2c2),
         ex6.diff(b2c1),
         ex6.diff(b2c2),
         ex6.diff(g2c1),
         ex6.diff(g2c2),
         ex6.diff(ad1c1),
         ex6.diff(ad1c2),
         ex6.diff(bd1c1),
         ex6.diff(bd1c2),
         ex6.diff(gd1c1),
         ex6.diff(gd1c2),
         ex6.diff(ad2c1),
         ex6.diff(ad2c2),
         ex6.diff(bd2c1),
         ex6.diff(bd2c2),
         ex6.diff(gd2c1),
         ex6.diff(gd2c2),
         ex6.diff(lamc1),
         ex6.diff(lamc2)],
        [ex7.diff(a1c1),
         ex7.diff(a1c2),
         ex7.diff(b1c1),
         ex7.diff(b1c2),
         ex7.diff(g1c1),
         ex7.diff(g1c2),
         ex7.diff(a2c1),
         ex7.diff(a2c2),
         ex7.diff(b2c1),
         ex7.diff(b2c2),
         ex7.diff(g2c1),
         ex7.diff(g2c2),
         ex7.diff(ad1c1),
         ex7.diff(ad1c2),
         ex7.diff(bd1c1),
         ex7.diff(bd1c2),
         ex7.diff(gd1c1),
         ex7.diff(gd1c2),
         ex7.diff(ad2c1),
         ex7.diff(ad2c2),
         ex7.diff(bd2c1),
         ex7.diff(bd2c2),
         ex7.diff(gd2c1),
         ex7.diff(gd2c2),
         ex7.diff(lamc1),
         ex7.diff(lamc2)],
        [ex8.diff(a1c1),
         ex8.diff(a1c2),
         ex8.diff(b1c1),
         ex8.diff(b1c2),
         ex8.diff(g1c1),
         ex8.diff(g1c2),
         ex8.diff(a2c1),
         ex8.diff(a2c2),
         ex8.diff(b2c1),
         ex8.diff(b2c2),
         ex8.diff(g2c1),
         ex8.diff(g2c2),
         ex8.diff(ad1c1),
         ex8.diff(ad1c2),
         ex8.diff(bd1c1),
         ex8.diff(bd1c2),
         ex8.diff(gd1c1),
         ex8.diff(gd1c2),
         ex8.diff(ad2c1),
         ex8.diff(ad2c2),
         ex8.diff(bd2c1),
         ex8.diff(bd2c2),
         ex8.diff(gd2c1),
         ex8.diff(gd2c2),
         ex8.diff(lamc1),
         ex8.diff(lamc2)],
        [ex9.diff(a1c1),
         ex9.diff(a1c2),
         ex9.diff(b1c1),
         ex9.diff(b1c2),
         ex9.diff(g1c1),
         ex9.diff(g1c2),
         ex9.diff(a2c1),
         ex9.diff(a2c2),
         ex9.diff(b2c1),
         ex9.diff(b2c2),
         ex9.diff(g2c1),
         ex9.diff(g2c2),
         ex9.diff(ad1c1),
         ex9.diff(ad1c2),
         ex9.diff(bd1c1),
         ex9.diff(bd1c2),
         ex9.diff(gd1c1),
         ex9.diff(gd1c2),
         ex9.diff(ad2c1),
         ex9.diff(ad2c2),
         ex9.diff(bd2c1),
         ex9.diff(bd2c2),
         ex9.diff(gd2c1),
         ex9.diff(gd2c2),
         ex9.diff(lamc1),
         ex9.diff(lamc2)],
        [ex10.diff(a1c1),
         ex10.diff(a1c2),
         ex10.diff(b1c1),
         ex10.diff(b1c2),
         ex10.diff(g1c1),
         ex10.diff(g1c2),
         ex10.diff(a2c1),
         ex10.diff(a2c2),
         ex10.diff(b2c1),
         ex10.diff(b2c2),
         ex10.diff(g2c1),
         ex10.diff(g2c2),
         ex10.diff(ad1c1),
         ex10.diff(ad1c2),
         ex10.diff(bd1c1),
         ex10.diff(bd1c2),
         ex10.diff(gd1c1),
         ex10.diff(gd1c2),
         ex10.diff(ad2c1),
         ex10.diff(ad2c2),
         ex10.diff(bd2c1),
         ex10.diff(bd2c2),
         ex10.diff(gd2c1),
         ex10.diff(gd2c2),
         ex10.diff(lamc1),
         ex10.diff(lamc2)],
        [ex11.diff(a1c1),
         ex11.diff(a1c2),
         ex11.diff(b1c1),
         ex11.diff(b1c2),
         ex11.diff(g1c1),
         ex11.diff(g1c2),
         ex11.diff(a2c1),
         ex11.diff(a2c2),
         ex11.diff(b2c1),
         ex11.diff(b2c2),
         ex11.diff(g2c1),
         ex11.diff(g2c2),
         ex11.diff(ad1c1),
         ex11.diff(ad1c2),
         ex11.diff(bd1c1),
         ex11.diff(bd1c2),
         ex11.diff(gd1c1),
         ex11.diff(gd1c2),
         ex11.diff(ad2c1),
         ex11.diff(ad2c2),
         ex11.diff(bd2c1),
         ex11.diff(bd2c2),
         ex11.diff(gd2c1),
         ex11.diff(gd2c2),
         ex11.diff(lamc1),
         ex11.diff(lamc2)],
        [ex12.diff(a1c1),
         ex12.diff(a1c2),
         ex12.diff(b1c1),
         ex12.diff(b1c2),
         ex12.diff(g1c1),
         ex12.diff(g1c2),
         ex12.diff(a2c1),
         ex12.diff(a2c2),
         ex12.diff(b2c1),
         ex12.diff(b2c2),
         ex12.diff(g2c1),
         ex12.diff(g2c2),
         ex12.diff(ad1c1),
         ex12.diff(ad1c2),
         ex12.diff(bd1c1),
         ex12.diff(bd1c2),
         ex12.diff(gd1c1),
         ex12.diff(gd1c2),
         ex12.diff(ad2c1),
         ex12.diff(ad2c2),
         ex12.diff(bd2c1),
         ex12.diff(bd2c2),
         ex12.diff(gd2c1),
         ex12.diff(gd2c2),
         ex12.diff(lamc1),
         ex12.diff(lamc2)],
        [ex13.diff(a1c1),
         ex13.diff(a1c2),
         ex13.diff(b1c1),
         ex13.diff(b1c2),
         ex13.diff(g1c1),
         ex13.diff(g1c2),
         ex13.diff(a2c1),
         ex13.diff(a2c2),
         ex13.diff(b2c1),
         ex13.diff(b2c2),
         ex13.diff(g2c1),
         ex13.diff(g2c2),
         ex13.diff(ad1c1),
         ex13.diff(ad1c2),
         ex13.diff(bd1c1),
         ex13.diff(bd1c2),
         ex13.diff(gd1c1),
         ex13.diff(gd1c2),
         ex13.diff(ad2c1),
         ex13.diff(ad2c2),
         ex13.diff(bd2c1),
         ex13.diff(bd2c2),
         ex13.diff(gd2c1),
         ex13.diff(gd2c2),
         ex13.diff(lamc1),
         ex13.diff(lamc2)],
        [ex14.diff(a1c1),
         ex14.diff(a1c2),
         ex14.diff(b1c1),
         ex14.diff(b1c2),
         ex14.diff(g1c1),
         ex14.diff(g1c2),
         ex14.diff(a2c1),
         ex14.diff(a2c2),
         ex14.diff(b2c1),
         ex14.diff(b2c2),
         ex14.diff(g2c1),
         ex14.diff(g2c2),
         ex14.diff(ad1c1),
         ex14.diff(ad1c2),
         ex14.diff(bd1c1),
         ex14.diff(bd1c2),
         ex14.diff(gd1c1),
         ex14.diff(gd1c2),
         ex14.diff(ad2c1),
         ex14.diff(ad2c2),
         ex14.diff(bd2c1),
         ex14.diff(bd2c2),
         ex14.diff(gd2c1),
         ex14.diff(gd2c2),
         ex14.diff(lamc1),
         ex14.diff(lamc2)],
        [ex15.diff(a1c1),
         ex15.diff(a1c2),
         ex15.diff(b1c1),
         ex15.diff(b1c2),
         ex15.diff(g1c1),
         ex15.diff(g1c2),
         ex15.diff(a2c1),
         ex15.diff(a2c2),
         ex15.diff(b2c1),
         ex15.diff(b2c2),
         ex15.diff(g2c1),
         ex15.diff(g2c2),
         ex15.diff(ad1c1),
         ex15.diff(ad1c2),
         ex15.diff(bd1c1),
         ex15.diff(bd1c2),
         ex15.diff(gd1c1),
         ex15.diff(gd1c2),
         ex15.diff(ad2c1),
         ex15.diff(ad2c2),
         ex15.diff(bd2c1),
         ex15.diff(bd2c2),
         ex15.diff(gd2c1),
         ex15.diff(gd2c2),
         ex15.diff(lamc1),
         ex15.diff(lamc2)],
        [ex16.diff(a1c1),
         ex16.diff(a1c2),
         ex16.diff(b1c1),
         ex16.diff(b1c2),
         ex16.diff(g1c1),
         ex16.diff(g1c2),
         ex16.diff(a2c1),
         ex16.diff(a2c2),
         ex16.diff(b2c1),
         ex16.diff(b2c2),
         ex16.diff(g2c1),
         ex16.diff(g2c2),
         ex16.diff(ad1c1),
         ex16.diff(ad1c2),
         ex16.diff(bd1c1),
         ex16.diff(bd1c2),
         ex16.diff(gd1c1),
         ex16.diff(gd1c2),
         ex16.diff(ad2c1),
         ex16.diff(ad2c2),
         ex16.diff(bd2c1),
         ex16.diff(bd2c2),
         ex16.diff(gd2c1),
         ex16.diff(gd2c2),
         ex16.diff(lamc1),
         ex16.diff(lamc2)],
        [ex17.diff(a1c1),
         ex17.diff(a1c2),
         ex17.diff(b1c1),
         ex17.diff(b1c2),
         ex17.diff(g1c1),
         ex17.diff(g1c2),
         ex17.diff(a2c1),
         ex17.diff(a2c2),
         ex17.diff(b2c1),
         ex17.diff(b2c2),
         ex17.diff(g2c1),
         ex17.diff(g2c2),
         ex17.diff(ad1c1),
         ex17.diff(ad1c2),
         ex17.diff(bd1c1),
         ex17.diff(bd1c2),
         ex17.diff(gd1c1),
         ex17.diff(gd1c2),
         ex17.diff(ad2c1),
         ex17.diff(ad2c2),
         ex17.diff(bd2c1),
         ex17.diff(bd2c2),
         ex17.diff(gd2c1),
         ex17.diff(gd2c2),
         ex17.diff(lamc1),
         ex17.diff(lamc2)],
        [ex18.diff(a1c1),
         ex18.diff(a1c2),
         ex18.diff(b1c1),
         ex18.diff(b1c2),
         ex18.diff(g1c1),
         ex18.diff(g1c2),
         ex18.diff(a2c1),
         ex18.diff(a2c2),
         ex18.diff(b2c1),
         ex18.diff(b2c2),
         ex18.diff(g2c1),
         ex18.diff(g2c2),
         ex18.diff(ad1c1),
         ex18.diff(ad1c2),
         ex18.diff(bd1c1),
         ex18.diff(bd1c2),
         ex18.diff(gd1c1),
         ex18.diff(gd1c2),
         ex18.diff(ad2c1),
         ex18.diff(ad2c2),
         ex18.diff(bd2c1),
         ex18.diff(bd2c2),
         ex18.diff(gd2c1),
         ex18.diff(gd2c2),
         ex18.diff(lamc1),
         ex18.diff(lamc2)],
        [ex19.diff(a1c1),
         ex19.diff(a1c2),
         ex19.diff(b1c1),
         ex19.diff(b1c2),
         ex19.diff(g1c1),
         ex19.diff(g1c2),
         ex19.diff(a2c1),
         ex19.diff(a2c2),
         ex19.diff(b2c1),
         ex19.diff(b2c2),
         ex19.diff(g2c1),
         ex19.diff(g2c2),
         ex19.diff(ad1c1),
         ex19.diff(ad1c2),
         ex19.diff(bd1c1),
         ex19.diff(bd1c2),
         ex19.diff(gd1c1),
         ex19.diff(gd1c2),
         ex19.diff(ad2c1),
         ex19.diff(ad2c2),
         ex19.diff(bd2c1),
         ex19.diff(bd2c2),
         ex19.diff(gd2c1),
         ex19.diff(gd2c2),
         ex19.diff(lamc1),
         ex19.diff(lamc2)],
        [ex20.diff(a1c1),
         ex20.diff(a1c2),
         ex20.diff(b1c1),
         ex20.diff(b1c2),
         ex20.diff(g1c1),
         ex20.diff(g1c2),
         ex20.diff(a2c1),
         ex20.diff(a2c2),
         ex20.diff(b2c1),
         ex20.diff(b2c2),
         ex20.diff(g2c1),
         ex20.diff(g2c2),
         ex20.diff(ad1c1),
         ex20.diff(ad1c2),
         ex20.diff(bd1c1),
         ex20.diff(bd1c2),
         ex20.diff(gd1c1),
         ex20.diff(gd1c2),
         ex20.diff(ad2c1),
         ex20.diff(ad2c2),
         ex20.diff(bd2c1),
         ex20.diff(bd2c2),
         ex20.diff(gd2c1),
         ex20.diff(gd2c2),
         ex20.diff(lamc1),
         ex20.diff(lamc2)],
        [ex21.diff(a1c1),
         ex21.diff(a1c2),
         ex21.diff(b1c1),
         ex21.diff(b1c2),
         ex21.diff(g1c1),
         ex21.diff(g1c2),
         ex21.diff(a2c1),
         ex21.diff(a2c2),
         ex21.diff(b2c1),
         ex21.diff(b2c2),
         ex21.diff(g2c1),
         ex21.diff(g2c2),
         ex21.diff(ad1c1),
         ex21.diff(ad1c2),
         ex21.diff(bd1c1),
         ex21.diff(bd1c2),
         ex21.diff(gd1c1),
         ex21.diff(gd1c2),
         ex21.diff(ad2c1),
         ex21.diff(ad2c2),
         ex21.diff(bd2c1),
         ex21.diff(bd2c2),
         ex21.diff(gd2c1),
         ex21.diff(gd2c2),
         ex21.diff(lamc1),
         ex21.diff(lamc2)],
        [ex22.diff(a1c1),
         ex22.diff(a1c2),
         ex22.diff(b1c1),
         ex22.diff(b1c2),
         ex22.diff(g1c1),
         ex22.diff(g1c2),
         ex22.diff(a2c1),
         ex22.diff(a2c2),
         ex22.diff(b2c1),
         ex22.diff(b2c2),
         ex22.diff(g2c1),
         ex22.diff(g2c2),
         ex22.diff(ad1c1),
         ex22.diff(ad1c2),
         ex22.diff(bd1c1),
         ex22.diff(bd1c2),
         ex22.diff(gd1c1),
         ex22.diff(gd1c2),
         ex22.diff(ad2c1),
         ex22.diff(ad2c2),
         ex22.diff(bd2c1),
         ex22.diff(bd2c2),
         ex22.diff(gd2c1),
         ex22.diff(gd2c2),
         ex22.diff(lamc1),
         ex22.diff(lamc2)],
        [ex23.diff(a1c1),
         ex23.diff(a1c2),
         ex23.diff(b1c1),
         ex23.diff(b1c2),
         ex23.diff(g1c1),
         ex23.diff(g1c2),
         ex23.diff(a2c1),
         ex23.diff(a2c2),
         ex23.diff(b2c1),
         ex23.diff(b2c2),
         ex23.diff(g2c1),
         ex23.diff(g2c2),
         ex23.diff(ad1c1),
         ex23.diff(ad1c2),
         ex23.diff(bd1c1),
         ex23.diff(bd1c2),
         ex23.diff(gd1c1),
         ex23.diff(gd1c2),
         ex23.diff(ad2c1),
         ex23.diff(ad2c2),
         ex23.diff(bd2c1),
         ex23.diff(bd2c2),
         ex23.diff(gd2c1),
         ex23.diff(gd2c2),
         ex23.diff(lamc1),
         ex23.diff(lamc2)],
        [ex24.diff(a1c1),
         ex24.diff(a1c2),
         ex24.diff(b1c1),
         ex24.diff(b1c2),
         ex24.diff(g1c1),
         ex24.diff(g1c2),
         ex24.diff(a2c1),
         ex24.diff(a2c2),
         ex24.diff(b2c1),
         ex24.diff(b2c2),
         ex24.diff(g2c1),
         ex24.diff(g2c2),
         ex24.diff(ad1c1),
         ex24.diff(ad1c2),
         ex24.diff(bd1c1),
         ex24.diff(bd1c2),
         ex24.diff(gd1c1),
         ex24.diff(gd1c2),
         ex24.diff(ad2c1),
         ex24.diff(ad2c2),
         ex24.diff(bd2c1),
         ex24.diff(bd2c2),
         ex24.diff(gd2c1),
         ex24.diff(gd2c2),
         ex24.diff(lamc1),
         ex24.diff(lamc2)],
        [ex25.diff(a1c1),
         ex25.diff(a1c2),
         ex25.diff(b1c1),
         ex25.diff(b1c2),
         ex25.diff(g1c1),
         ex25.diff(g1c2),
         ex25.diff(a2c1),
         ex25.diff(a2c2),
         ex25.diff(b2c1),
         ex25.diff(b2c2),
         ex25.diff(g2c1),
         ex25.diff(g2c2),
         ex25.diff(ad1c1),
         ex25.diff(ad1c2),
         ex25.diff(bd1c1),
         ex25.diff(bd1c2),
         ex25.diff(gd1c1),
         ex25.diff(gd1c2),
         ex25.diff(ad2c1),
         ex25.diff(ad2c2),
         ex25.diff(bd2c1),
         ex25.diff(bd2c2),
         ex25.diff(gd2c1),
         ex25.diff(gd2c2),
         ex25.diff(lamc1),
         ex25.diff(lamc2)],
        [ex26.diff(a1c1),
         ex26.diff(a1c2),
         ex26.diff(b1c1),
         ex26.diff(b1c2),
         ex26.diff(g1c1),
         ex26.diff(g1c2),
         ex26.diff(a2c1),
         ex26.diff(a2c2),
         ex26.diff(b2c1),
         ex26.diff(b2c2),
         ex26.diff(g2c1),
         ex26.diff(g2c2),
         ex26.diff(ad1c1),
         ex26.diff(ad1c2),
         ex26.diff(bd1c1),
         ex26.diff(bd1c2),
         ex26.diff(gd1c1),
         ex26.diff(gd1c2),
         ex26.diff(ad2c1),
         ex26.diff(ad2c2),
         ex26.diff(bd2c1),
         ex26.diff(bd2c2),
         ex26.diff(gd2c1),
         ex26.diff(gd2c2),
         ex26.diff(lamc1),
         ex26.diff(lamc2)]
        
    ])

conlist = np.array([
    ex1,
    ex2,
    ex3,
    ex4,
    ex5,
    ex6,
    ex7,
    ex8,
    ex9,
    ex10,
    ex11,
    ex12,
    ex13,
    ex14,
    ex15,
    ex16,
    ex17,
    ex18,
    ex19,
    ex20,
    ex21,
    ex22,
    ex23,
    ex24,
    ex25,
    ex26
])

# Solver Setup

# Time array based on time step
dt = .1    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Simulation data container

# H3
rs2simdatalist = np.zeros((t_arr.shape[0],12+1))

# Initial Data
v = 1.      # Initial Velocity
x = 1.      # Rod Length 
m = 1.      # Mass of point masses
params = [v,x,m]


# Sim bar in H3
startvec = np.array([
    .5,np.pi/2.,np.pi/2.,
    .5,np.pi/2.,3.*np.pi/2.,
    killingvech3([.5,np.pi/2.,np.pi/2.],v,"x")[0],
    killingvech3([.5,np.pi/2.,np.pi/2.],v,"x")[1],
    killingvech3([.5,np.pi/2.,np.pi/2.],v,"x")[2], 
    killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")[0],
    killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")[1],
    killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")[2],
    0.]).flatten()

# a1c1,b1c1,g1c1 = startvec[0]
# a1c2,b1c2,g1c2 = startvec[0]
# a2c1,b2c1,g2c1 = startvec[1]
# a2c2,b2c2,g2c2 = startvec[1]

# ad1c1,bd1c1,gd1c1 = startvec[2]
# ad1c2,bd1c2,gd1c2 = startvec[2]
# ad2c1,bd2c1,gd2c1 = startvec[3]
# ad2c2,bd2c2,gd2c2 = startvec[3]

# lamc1 = startvec[4]
# lamc2 = startvec[4]

# Sim bar in S3
# startvec = np.array([
#     [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
#     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()

# First Step
step = 0

# H3
rs2simdatalist[step] = startvec

initvec = [
     startvec[0],startvec[0],
     startvec[1],startvec[1],
     startvec[2],startvec[2],
     startvec[3],startvec[3],
     startvec[4],startvec[4],
     startvec[5],startvec[5],

     startvec[6],startvec[6],
     startvec[7],startvec[7],
     startvec[8],startvec[8],
     startvec[9],startvec[9],
     startvec[10],startvec[10],
     startvec[11],startvec[11],

     startvec[12],startvec[12]
     ]

# Update input
inputvals[a1n] = startvec[0]
inputvals[a1c1] = startvec[0]
inputvals[a1c2] = startvec[0]
inputvals[b1n] = startvec[1]
inputvals[b1c1] = startvec[1]
inputvals[b1c2] = startvec[1]
inputvals[g1n] = startvec[2]
inputvals[g1c1] = startvec[2]
inputvals[g1c2] = startvec[2]
inputvals[a2n] = startvec[3]
inputvals[a2c1] = startvec[3]
inputvals[a2c2] = startvec[3]
inputvals[b2n] = startvec[4]
inputvals[b2c1] = startvec[4]
inputvals[b2c2] = startvec[4]
inputvals[g2n] = startvec[5]
inputvals[g2c1] = startvec[5]
inputvals[g2c2] = startvec[5]

inputvals[ad1n] = startvec[6]
inputvals[ad1c1] = startvec[6]
inputvals[ad1c2] = startvec[6]
inputvals[bd1n] = startvec[7]
inputvals[bd1c1] = startvec[7]
inputvals[bd1c2] = startvec[7]
inputvals[gd1n] = startvec[8]
inputvals[gd1c1] = startvec[8]
inputvals[gd1c2] = startvec[8]
inputvals[ad2n] = startvec[9]
inputvals[ad2c1] = startvec[9]
inputvals[ad2c2] = startvec[9]
inputvals[bd2n] = startvec[10]
inputvals[bd2c1] = startvec[10]
inputvals[bd2c2] = startvec[10]
inputvals[gd2n] = startvec[11]
inputvals[gd2c1] = startvec[11]
inputvals[gd2c2] = startvec[11]

inputvals[lamc1] = startvec[12]
inputvals[lamc2] = startvec[12]


jmat1 = evaljac(jacobian,inputvals)
con1 = evalcon(conlist,inputvals)

diff1 = np.linalg.solve(jmat1.astype('float64'),-con1.astype('float64'))
print(np.linalg.norm(diff1))

val1 = diff1 + initvec

counter = 0
while (np.linalg.norm(diff1) >= 1e-8 and counter <= 100):
    startvec = val1[::2]
    # Update input
    inputvals[a1n] = startvec[0]
    inputvals[a1c1] = startvec[0]
    inputvals[a1c2] = startvec[0]
    inputvals[b1n] = startvec[1]
    inputvals[b1c1] = startvec[1]
    inputvals[b1c2] = startvec[1]
    inputvals[g1n] = startvec[2]
    inputvals[g1c1] = startvec[2]
    inputvals[g1c2] = startvec[2]
    inputvals[a2n] = startvec[3]
    inputvals[a2c1] = startvec[3]
    inputvals[a2c2] = startvec[3]
    inputvals[b2n] = startvec[4]
    inputvals[b2c1] = startvec[4]
    inputvals[b2c2] = startvec[4]
    inputvals[g2n] = startvec[5]
    inputvals[g2c1] = startvec[5]
    inputvals[g2c2] = startvec[5]

    inputvals[ad1n] = startvec[6]
    inputvals[ad1c1] = startvec[6]
    inputvals[ad1c2] = startvec[6]
    inputvals[bd1n] = startvec[7]
    inputvals[bd1c1] = startvec[7]
    inputvals[bd1c2] = startvec[7]
    inputvals[gd1n] = startvec[8]
    inputvals[gd1c1] = startvec[8]
    inputvals[gd1c2] = startvec[8]
    inputvals[ad2n] = startvec[9]
    inputvals[ad2c1] = startvec[9]
    inputvals[ad2c2] = startvec[9]
    inputvals[bd2n] = startvec[10]
    inputvals[bd2c1] = startvec[10]
    inputvals[bd2c2] = startvec[10]
    inputvals[gd2n] = startvec[11]
    inputvals[gd2c1] = startvec[11]
    inputvals[gd2c2] = startvec[11]

    inputvals[lamc1] = startvec[12]
    inputvals[lamc2] = startvec[12]

    # jmat2 = evaljac(jacobian,inputvals)
    con2 = evalcon(conlist,inputvals)

    diff2 = np.linalg.solve(jmat1.astype('float64'),-con2.astype('float64'))

    val2 = diff2 + val1

    val1 = val2
    diff1 = diff2
    print(np.linalg.norm(diff1))

    counter += 1
     

# Exact in H3
# gs1exactdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)


# Exact in H3
# startvecgs1ex = gs1exactdatalist[step]
# startvecgs2ex = gs2exactdatalist[step]
# startvecgs3ex = gs3exactdatalist[step]

# Sim in H3
# startvecgs1sim = gs1simdatalist[step]
# startvecgs2sim = gs2simdatalist[step]
# startvecgs3sim = gs3simdatalist[step]

# Exact in S3
# startvecgs1ex = gs1exactdatalist[step]
# startvecgs2ex = gs2exactdatalist[step]
# startvecgs3ex = gs3exactdatalist[step]

# Sim in S3
# startvecgs1sim = gs1simdatalist[step]
# startvecgs2sim = gs2simdatalist[step]
# startvecgs3sim = gs3simdatalist[step]

# step += 1

# while (step <= int(t_max/dt)):
#     # Exact in H3
#     # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)
#     # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)
#     # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)

#     # Sim in H3
#     # gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt)
#     # gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt)
#     # gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt)    

#     # Exact in S3
#     # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)
#     # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)
#     # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)

#     # Sim in S3
#     gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt)
#     gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt)
#     gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt) 

#     # Exact in H3
#     # startvecgs1ex = gs1exactdatalist[step]
#     # startvecgs2ex = gs2exactdatalist[step]
#     # startvecgs3ex = gs3exactdatalist[step]

#     # Sim in H3
#     # startvecgs1sim = gs1simdatalist[step]
#     # startvecgs2sim = gs2simdatalist[step]
#     # startvecgs3sim = gs3simdatalist[step]

#     # Exact in S3
#     # startvecgs1ex = gs1exactdatalist[step]
#     # startvecgs2ex = gs2exactdatalist[step]
#     # startvecgs3ex = gs3exactdatalist[step]

#     # Sim in S3
#     startvecgs1sim = gs1simdatalist[step]
#     startvecgs2sim = gs2simdatalist[step]
#     startvecgs3sim = gs3simdatalist[step]

#     if step%int(1/dt)==0:
#             print(step)
#     step += 1

# # Exact in H3
# # np.save("gausss1_tmax10_dt01",gs1exactdatalist)
# # np.save("gausss2_tmax10_dt01",gs2exactdatalist)
# # np.save("h3_r_gausss3_gt_tmax10_dt000005",gs3exactdatalist)

# # Sim in H3
# # np.save("s3_r_gausss1_sim_tmax10_dt1",gs1simdatalist)
# # np.save("s3_r_gausss2_sim_tmax10_dt1",gs2simdatalist)
# # np.save("s3_r_gausss3_sim_tmax10_dt1",gs3simdatalist)

# # Exact in S3
# # np.save("gausss1_tmax10_dt01",gs1exactdatalist)
# # np.save("gausss2_tmax10_dt01",gs2exactdatalist)
# # np.save("gausss3_gt_tmax10_dt000005",gs3exactdatalist)

# # Sim in S3
# np.save("s3_r_gausss1_sim_tmax10_dt0001",gs1simdatalist)
# np.save("s3_r_gausss2_sim_tmax10_dt0001",gs2simdatalist)
# np.save("s3_r_gausss3_sim_tmax10_dt0001",gs3simdatalist)

# data1 = np.load("s3_r_gausss1_sim_tmax10_dt0001.npy")
# data2 = np.load("s3_r_gausss2_sim_tmax10_dt0001.npy")
# data3 = np.load("s3_r_gausss3_sim_tmax10_dt0001.npy")

# # data1 = np.load("h3_r_gausss3_gt_tmax10_dt000005.npy")


# fig,ax=plt.subplots(1,1)

# distdata1 = np.zeros(np.shape(data1)[0])
# distdata2 = np.zeros(np.shape(data2)[0])
# distdata3 = np.zeros(np.shape(data3)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
#     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
#     distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
#     counter += 1

# # counter = 0
# # for a in range(np.shape(data1)[0]):
# #     distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
# #     distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
# #     distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
# #     counter += 1

# # ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# # ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# # ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")
# # ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
# ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()

