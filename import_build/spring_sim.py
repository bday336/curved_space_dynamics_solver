# Allow for package importing
import os, sys
sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.path.dirname(os.getcwd())+"/src")

import numpy as np
import matplotlib.pyplot as plt
from random import random,uniform

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import animation, rc

# Import Packages
from src.utils.function_bank import boostxh3, rotxh3, rotzh3, hyp2poin3d, rot2hyp

from src.Computation.DataList import DataList

from src.Computation.State import State

from src.ConfigurationSpace.ConfigurationSpace import ConfigurationSpace
from src.ConfigurationSpace.Simulation import Simulation

from src.AmbientSpace.ExampleSpaces.Euclidean import euclidean
from src.AmbientSpace.ExampleSpaces.HypSpaceCoords import hyperbolic
# from src.AmbientSpace.ExampleSpaces.SphericalStereoProj import spherical
# from src.AmbientSpace.ExampleSpaces.InhomogeneousNeg import inhomogeneousNeg
# from src.AmbientSpace.ExampleSpaces.InhomogeneousPos import inhomogeneousPos
# from src.AmbientSpace.ExampleSpaces.H2xE import h2xe
# from src.AmbientSpace.ExampleSpaces.S2xE import s2xe

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

# Set the ambient space for the simulation environment
ambientSpace = hyperbolic

# Build a configuration space
NumBalls = 2
MaxRad = .2 #ambientSpace.obstacle.size/5.

radii = []
masses = []
for a in range(NumBalls):
    r = MaxRad #* random()+0.05
    m = 1. #10.*r*r*r
    radii.append(r)
    masses.append(m)

stiffness = 1.
eq_len = 1.


configurationSpace = ConfigurationSpace(masses, radii, ambientSpace)
# // let maxPos = 2.;
# // let maxVel = 1;

# //build the initial set of states for the system:

# E3
# iniCond = [
#     State(np.array([4,0,0]),np.array([-1,0,0])),
#     State(np.array([0,4,.1]),np.array([0,-1,0]))
# ]

# H3
iniCond = [
    State(np.array([.5,np.pi/2.,0]),np.array([0,0,0])),
    State(np.array([.5,np.pi/2.,np.pi]),np.array([0,0,-1]))
]

dataList = DataList(iniCond)

dt = 0.001
tmaxNum = 100000

# //make the simulation
sim = Simulation( ambientSpace, dataList, configurationSpace, dt )

for c in range(tmaxNum):
    sim.step()

# //make the visualization of the simulation
# let viz = new RenderSim( sim, radii );
    
## Plot Space Trajectory and Distance Error
fig = plt.figure(figsize=(6,6))

## Plot the Space Trajectory of the Rod Body System
ax1=fig.add_subplot(1,1,1,projection='3d')

# E3

# # Draw Unit Sphere Horizon in Embedding Space (using Poincare Disk Model of Hyperbolic Space)
# u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
# x = 6.*np.sin(u)*np.cos(v)
# y = 6.*np.sin(u)*np.sin(v)
# z = 6.*np.cos(u)
# ax1.plot_wireframe(x, y, z, color="b", alpha=.1)
# ax1.set_xlim3d(-10,10)
# ax1.set_xlabel('X')
# ax1.set_ylim3d(-10,10)
# ax1.set_ylabel('Y')
# ax1.set_zlim3d(-10,10)
# ax1.set_zlabel('Z')

# H3

# Draw Unit Sphere Horizon in Embedding Space (using Poincare Disk Model of Hyperbolic Space)
u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
x = 1.*np.sin(u)*np.cos(v)
y = 1.*np.sin(u)*np.sin(v)
z = 1.*np.cos(u)
ax1.plot_wireframe(x, y, z, color="b", alpha=.1)
ax1.set_xlim3d(-1,1)
ax1.set_xlabel('X')
ax1.set_ylim3d(-1,1)
ax1.set_ylabel('Y')
ax1.set_zlim3d(-1,1)
ax1.set_zlabel('Z')

ub, vb = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
xb = np.sinh(2)/(np.cosh(2)+1) * np.sin(ub)*np.cos(vb)
yb = np.sinh(2)/(np.cosh(2)+1) * np.sin(ub)*np.sin(vb)
zb = np.sinh(2)/(np.cosh(2)+1) * np.cos(ub)
ax1.plot_wireframe(xb, yb, zb, color="k", alpha=.1)

hyppart1plot=[]
hyppart2plot=[]
part1plot=[]
part2plot=[]

# E3
# for a in sim.data_container:
#     part1plot.append(a.data[0].pos)
#     part2plot.append(a.data[1].pos)

# part1plot = np.array(part1plot)
# part2plot = np.array(part2plot)
    
# # # Generate the pseudo-shell surface of the vertices of the rod body
# part1x,part1y,part1z=genballe3(part1plot[-1],radii[0])
# part2x,part2y,part2z=genballe3(part2plot[-1],radii[1])

#H3
for a in sim.data_container:
    hyppart1plot.append(list(rot2hyp(a.data[0].pos)))
    hyppart2plot.append(list(rot2hyp(a.data[1].pos)))

    part1plot.append(list(hyp2poin3d(rot2hyp(a.data[0].pos))))
    part2plot.append(list(hyp2poin3d(rot2hyp(a.data[1].pos))))

hyppart1plot = np.array(hyppart1plot)
hyppart2plot = np.array(hyppart2plot)
part1plot = np.array(part1plot)
part2plot = np.array(part2plot)

# Generate the pseudo-shell surface of the vertices of the rod body
part1x,part1y,part1z=genballh3(hyppart1plot[-1],radii[0])
part2x,part2y,part2z=genballh3(hyppart2plot[-1],radii[1])


# # Draw the pseudo-shell surface of the vertices of the rod body
ax1.plot_surface(part1x, part1y, part1z, color="b")
ax1.plot_surface(part2x, part2y, part2z, color="b")


# Draw trajectory of vertices of rod body
ax1.plot3D(part1plot[:,0],part1plot[:,1],part1plot[:,2], label="particle 1")
ax1.plot3D(part2plot[:,0],part2plot[:,1],part2plot[:,2], label="particle 2")
ax1.legend(loc= 'lower left')

# ## Plot the error in distance between vertices compared to fixed distance
# ax2=fig.add_subplot(1,2,2)

# # Generate distance error data
# distdata1 = np.zeros(np.shape(data1)[0])

# counter = 0
# for b in range(np.shape(data1)[0]):
#     distdata1[counter] = h3dist(rot2hyp(data1[b][0:3]),rot2hyp(data1[b][3:6]))
#     counter += 1

# # Draw distance error as a function of simulation time
# ax2.plot(sim_test.t_arr,distdata1,marker = ".",color='k',linestyle = "None",label = "Gauss s1")
# ax2.plot(sim_test.t_arr,np.full(sim_test.t_arr.shape,(r1+r2)),color='r',label = "Collision Boundary")
# ax2.legend()
# ax2.set_title('Simulation Data')
# ax2.set_xlabel('t')
# ax2.set_ylabel('l')

fig.tight_layout()	

plt.show()

# ------------------------------------------------------------------
### Uncomment to just generate gif of trajectory of the particle ###
# ------------------------------------------------------------------

# # Generate gif
# # create empty lists for the x and y data
# x1 = []
# y1 = []
# z1 = []
# x2 = []
# y2 = []
# z2 = []
# x3 = []
# y3 = []
# z3 = []

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(111, projection='3d')
# ax1.set_aspect("equal")

# #draw sphere
# Draw Unit Sphere Horizon in Embedding Space (using Poincare Disk Model of Hyperbolic Space)
u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
x = 1.*np.sin(u)*np.cos(v)
y = 1.*np.sin(u)*np.sin(v)
z = 1.*np.cos(u)
ax1.plot_wireframe(x, y, z, color="b", alpha=.1)
ax1.set_xlim3d(-1,1)
ax1.set_xlabel('X')
ax1.set_ylim3d(-1,1)
ax1.set_ylabel('Y')
ax1.set_zlim3d(-1,1)
ax1.set_zlabel('Z')

ub, vb = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
xb = np.sinh(2)/(np.cosh(2)+1) * np.sin(ub)*np.cos(vb)
yb = np.sinh(2)/(np.cosh(2)+1) * np.sin(ub)*np.sin(vb)
zb = np.sinh(2)/(np.cosh(2)+1) * np.cos(ub)
ax1.plot_wireframe(xb, yb, zb, color="k", alpha=.1)

# Particle Plot data
ball_box=[]
viz_balls=[]

ball_box.append(genballh3(hyppart1plot[0],radii[0])) #hypercirch3(array([sinh(gat[a])*sin(gbt[a])*cos(ggt[a]),sinh(gat[a])*sin(gbt[a])*sin(ggt[a]),sinh(gat[a])*cos(gbt[a]),cosh(gat[a])]),.1))
viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")])

ball_box.append(genballh3(hyppart2plot[0],radii[0])) #hypercirch3(array([sinh(gat[a])*sin(gbt[a])*cos(ggt[a]),sinh(gat[a])*sin(gbt[a])*sin(ggt[a]),sinh(gat[a])*cos(gbt[a]),cosh(gat[a])]),.1))
viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")])

# #draw trajectory
# for b in range(vert_num):
#     ax1.plot3D(gut[b::vert_num],gvt[b::vert_num],grt[b::vert_num], label="particle []".format(b))
# ax1.legend(loc= 'lower left')

# part1x,part1y,part1z=hypercirch3(array([sinh(gat[0::3][0])*sin(gbt[0::3][0])*cos(ggt[0::3][0]),sinh(gat[0::3][0])*sin(gbt[0::3][0])*sin(ggt[0::3][0]),sinh(gat[0::3][0])*cos(gbt[0::3][0]),cosh(gat[0::3][0])]),particles[0][7])
# part2x,part2y,part2z=hypercirch3(array([sinh(gat[1::3][1])*sin(gbt[1::3][1])*cos(ggt[1::3][1]),sinh(gat[1::3][1])*sin(gbt[1::3][1])*sin(ggt[1::3][1]),sinh(gat[1::3][1])*cos(gbt[1::3][1]),cosh(gat[1::3][1])]),particles[1][7])
# part3x,part3y,part3z=hypercirch3(array([sinh(gat[2::3][2])*sin(gbt[2::3][2])*cos(ggt[2::3][2]),sinh(gat[2::3][2])*sin(gbt[2::3][2])*sin(ggt[2::3][2]),sinh(gat[2::3][2])*cos(gbt[2::3][2]),cosh(gat[2::3][2])]),particles[2][7])
# ball1=[ax1.plot_surface(part1x,part1y,part1z, color="b")]
# ball2=[ax1.plot_surface(part2x,part2y,part2z, color="r")]
# ball3=[ax1.plot_surface(part3x,part3y,part3z, color="k")]

# animation function. This is called sequentially
frames=100
def animate(i):
    ax1.plot3D(part1plot[:int(tmaxNum*i/frames),0],part1plot[:int(tmaxNum*i/frames),1],part1plot[:int(tmaxNum*i/frames),2],)
    ax1.plot3D(part2plot[:int(tmaxNum*i/frames),0],part2plot[:int(tmaxNum*i/frames),1],part2plot[:int(tmaxNum*i/frames),2],)
    # ax1.legend(loc= 'lower left')

    ball_box.append(genballh3(hyppart1plot[int(tmaxNum*i/frames)],radii[0])) #array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
    viz_balls[0][0].remove()
    viz_balls[0][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")
    ball_box.append(genballh3(hyppart2plot[int(tmaxNum*i/frames)],radii[1])) #array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
    viz_balls[1][0].remove()
    viz_balls[1][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")
    
    # ax1.plot3D(gut[b::vert_num][:int(len(timearr)*i/frames)],gvt[b::vert_num][:int(len(timearr)*i/frames)],grt[b::vert_num][:int(len(timearr)*i/frames)], label="particle 1")
    # for b in range(vert_num):
    #     ax1.plot3D(gut[b::vert_num][:int(len(timearr)*i/frames)],gvt[b::vert_num][:int(len(timearr)*i/frames)],grt[b::vert_num][:int(len(timearr)*i/frames)], label="particle {}".format(b))
    # ax1.plot3D(gut[0::3][:int(len(timearr)*i/frames)],gvt[0::3][:int(len(timearr)*i/frames)],grt[0::3][:int(len(timearr)*i/frames)])
    # ax1.plot3D(gut[1::3][:int(len(timearr)*i/frames)],gvt[1::3][:int(len(timearr)*i/frames)],grt[1::3][:int(len(timearr)*i/frames)])
    # ax1.plot3D(gut[2::3][:int(len(timearr)*i/frames)],gvt[2::3][:int(len(timearr)*i/frames)],grt[2::3][:int(len(timearr)*i/frames)])
        # part_plot.append(genballh3(array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
        # viz_balls[b][0].remove()
        # viz_balls[b][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")
    # part1x,part1y,part1z=hypercirch3(array([sinh(gat[0::3][int(len(timearr)*i/frames)])*sin(gbt[0::3][int(len(timearr)*i/frames)])*cos(ggt[0::3][int(len(timearr)*i/frames)]),sinh(gat[0::3][int(len(timearr)*i/frames)])*sin(gbt[0::3][int(len(timearr)*i/frames)])*sin(ggt[0::3][int(len(timearr)*i/frames)]),sinh(gat[0::3][int(len(timearr)*i/frames)])*cos(gbt[0::3][int(len(timearr)*i/frames)]),cosh(gat[0::3][int(len(timearr)*i/frames)])]),particles[0][7])
    # part2x,part2y,part2z=hypercirch3(array([sinh(gat[1::3][int(len(timearr)*i/frames)])*sin(gbt[1::3][int(len(timearr)*i/frames)])*cos(ggt[1::3][int(len(timearr)*i/frames)]),sinh(gat[1::3][int(len(timearr)*i/frames)])*sin(gbt[1::3][int(len(timearr)*i/frames)])*sin(ggt[1::3][int(len(timearr)*i/frames)]),sinh(gat[1::3][int(len(timearr)*i/frames)])*cos(gbt[1::3][int(len(timearr)*i/frames)]),cosh(gat[1::3][int(len(timearr)*i/frames)])]),particles[1][7])
    # part3x,part3y,part3z=hypercirch3(array([sinh(gat[2::3][int(len(timearr)*i/frames)])*sin(gbt[2::3][int(len(timearr)*i/frames)])*cos(ggt[2::3][int(len(timearr)*i/frames)]),sinh(gat[2::3][int(len(timearr)*i/frames)])*sin(gbt[2::3][int(len(timearr)*i/frames)])*sin(ggt[2::3][int(len(timearr)*i/frames)]),sinh(gat[2::3][int(len(timearr)*i/frames)])*cos(gbt[2::3][int(len(timearr)*i/frames)]),cosh(gat[2::3][int(len(timearr)*i/frames)])]),particles[2][7])
    # ball1[0].remove()
    # ball1[0]=ax1.plot_surface(part1x,part1y,part1z, color="b")
    # ball2[0].remove()
    # ball2[0]=ax1.plot_surface(part2x,part2y,part2z, color="r")
    # ball3[0].remove()
    # ball3[0]=ax1.plot_surface(part3x,part3y,part3z, color="k")

# equivalent to rcParams['animation.html'] = 'html5'
rc('animation', html='html5')

# call the animator. blit=True means only re-draw the parts that 
# have changed.
anim = animation.FuncAnimation(fig, animate,frames=frames, interval=50)

anim.save('./h3springtest.mov', writer='imagemagick')