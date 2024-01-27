from sympy import *
import sympy as sp
import numpy as np

t,lamz = sp.symbols('t lamz')
a1 = sp.Function('a1')(t)
a2 = sp.Function('a2')(t)

ad1 = a1.diff(t)
ad2 = a2.diff(t)

lamx = sp.Function('lamx')(t)
lamy = sp.Function('lamy')(t)

lamxd = lamx.diff(t)
lamyd = lamy.diff(t)

# ad1 = sp.Function('ad1')(t)
# ad2 = sp.Function('ad2')(t)
# a1,a2 = sp.symbols('a1 a2')
# ad1,ad2 = sp.symbols('ad1 ad2')

kt,l,dt = sp.symbols('kt,l,dt')

def inadj(x,y,th):
    return Matrix([
        [cos(th) ,sin(th),x*sin(th)-y*cos(th)],
        [-sin(th),cos(th),x*cos(th)+y*sin(th)],
        [0,0,1]
        ])

je1 = Matrix([
    [cos(a1),-sin(a1),l/2*sin(a1)     ,0  ,0],
    [sin(a1),cos(a1) ,-l/2*(1+cos(a1)),l/2,0],
    [0,0,1,-1,0]
])

je2 = Matrix([
    [1,0,0,0,0],
    [0,1,0,0,0],
    [0,0,1,0,0]
])

je3 = Matrix([
    [cos(a2) ,sin(a2),l/2*sin(a2)    ,0,0],
    [-sin(a2),cos(a2),l/2*(1+cos(a2)),0,l/2],
    [0,0,1,0,1]
])

dadj1 = (inadj(-l/2,0,0) @ inadj(0,0,-a1) @ inadj(-l/2,0,0)).transpose()
dadj3 = (inadj(l/2,0,0) @ inadj(0,0,a2) @ inadj(l/2,0,0)).transpose()

dmat = Matrix([
    [kt*l,0     ,0         ],
    [0   ,2*kt*l,0         ],
    [0   ,0.    ,kt/12*l**3]
    ])

omg = dadj1 @ dmat @ je1 + dmat @ je2 + dadj3 @ dmat @ je3

acon = simplify(omg[0:3,0:3].inv() @ omg[0:3,3:6])


# Using the right trivialization of the lie algrebra
curvbx = Matrix([
    [
     acon[0,0].diff(a1) - acon[0,0].diff(a1) - (acon[1,0]*acon[2,0] - acon[2,0]*acon[1,0]),
     acon[0,0].diff(a2) - acon[0,1].diff(a1) - (acon[1,0]*acon[2,1] - acon[2,0]*acon[1,1])
    ],
    [
     acon[0,1].diff(a1) - acon[0,0].diff(a2) - (acon[1,1]*acon[2,0] - acon[2,1]*acon[1,0]),
     acon[0,1].diff(a2) - acon[0,1].diff(a2) - (acon[1,1]*acon[2,1] - acon[2,1]*acon[1,1])
    ]
])

curvby = Matrix([
    [
     acon[1,0].diff(a1) - acon[1,0].diff(a1) - (-acon[0,0]*acon[2,0] + acon[2,0]*acon[0,0]),
     acon[1,0].diff(a2) - acon[1,1].diff(a1) - (-acon[0,0]*acon[2,1] + acon[2,0]*acon[0,1])
    ],
    [
     acon[1,1].diff(a1) - acon[1,0].diff(a2) - (-acon[0,1]*acon[2,0] + acon[2,1]*acon[0,0]),
     acon[1,1].diff(a2) - acon[1,1].diff(a2) - (-acon[0,1]*acon[2,1] + acon[2,1]*acon[0,1])
    ]
])

curvbth = Matrix([
    [
     acon[2,0].diff(a1) - acon[2,0].diff(a1),
     acon[2,0].diff(a2) - acon[2,1].diff(a1)
    ],
    [
     acon[2,1].diff(a1) - acon[2,0].diff(a2),
     acon[2,1].diff(a2) - acon[2,1].diff(a2)
    ]
])

j1a = je1 @ Matrix([
    [-acon[0,0],-acon[0,1]],
    [-acon[1,0],-acon[1,1]],
    [-acon[2,0],-acon[2,1]],
    [1         ,0],
    [0         ,1]
])

j2a = je2 @ Matrix([
    [-acon[0,0],-acon[0,1]],
    [-acon[1,0],-acon[1,1]],
    [-acon[2,0],-acon[2,1]],
    [1         ,0],
    [0         ,1]
])

j3a = je3 @ Matrix([
    [-acon[0,0],-acon[0,1]],
    [-acon[1,0],-acon[1,1]],
    [-acon[2,0],-acon[2,1]],
    [1         ,0],
    [0         ,1]
])

d1a = (j1a.transpose()) @ dmat @ j1a
d2a = (j2a.transpose()) @ dmat @ j2a
d3a = (j3a.transpose()) @ dmat @ j3a

pterm = Matrix([ [ad1,ad2] ]) @ (d1a + d2a + d3a) @ Matrix([ [ad1], [ad2] ])


eqn1 = ((pterm[0,0].diff(ad1)).diff(t)) - pterm[0,0].diff(a1) + lamx*(curvbx[0, 0]*ad1 + curvbx[0, 1]*ad2) + lamy*(curvby[0, 0]*ad1 + curvby[0, 1]*ad2) + lamz*(curvbth[0, 0]*ad1 + curvbth[0, 1]*ad2)
eqn2 = ((pterm[0,0].diff(ad2)).diff(t)) - pterm[0,0].diff(a2) + lamx*(curvbx[1, 0]*ad1 + curvbx[1, 1]*ad2) + lamy*(curvby[1, 0]*ad1 + curvby[1, 1]*ad2) + lamz*(curvbth[1, 0]*ad1 + curvbth[1, 1]*ad2)

eqn3 = lamx.diff(t) + lamy*(acon[2, 0]*ad1 + acon[2, 1]*ad2)
eqn4 = lamy.diff(t) - lamx*(acon[2, 0]*ad1 + acon[2, 1]*ad2)

print("solving the system...")

sol = solve([eqn1,eqn2,eqn3,eqn4], ad1.diff(t), ad2.diff(t), lamx.diff(t), lamy.diff(t), dict=True)

