# Allow for package importing
import os, sys
sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.path.dirname(os.getcwd())+"/src")

import numpy as np
import matplotlib.pyplot as plt
from random import random,uniform

np.set_printoptions(linewidth=260) 
np.set_printoptions(precision=1)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import animation, rc

# Import Packages
from src.utils.mesh_bank import particle_mesh, dumbbell_mesh

from src.utils.function_bank import boostxh3, rotxh3, rotzh3, hyp2poin3d, rot2hyp, genballe3, genballh3

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

# Set the ambient space for the simulation environment
ambientSpace = hyperbolic

# Generate System
# //build the initial set of states for the system:

# E3
# iniCond = [
#     State(np.array([4,0,0]),np.array([-1,0,0])),
#     State(np.array([0,4,.1]),np.array([0,-1,0]))
# ]

# H3
# iniCond = [
#     State(np.array([.5,np.pi/2.,0]),np.array([0,0,0])),
#     State(np.array([.5,np.pi/2.,np.pi]),np.array([0,0,-1]))
# ]

# sim_system = dumbbell_mesh(
#     [np.array([.5,np.pi/2.,0*np.pi/2.]),np.array([.5,np.pi/2.,2.*np.pi/2.])],      # Positions
#     [np.array([0,0,-2.16395]),np.array([0,0,2.16395])],                                         # Velocities
#     isRigid=True)

# sim_system = dumbbell_mesh(
#     [np.array([.5,np.pi/2.,1*np.pi/2.]),np.array([.5,np.pi/2.,3.*np.pi/2.])],      # Positions
#     [np.array([0,0,-2.16395]),np.array([0,0,2.16395])],                                         # Velocities
#     isRigid=True).combine(
#         dumbbell_mesh(
#             [np.array([.5,np.pi/2.,0.*np.pi/2.]),np.array([.5,np.pi/2.,2.*np.pi/2.])],      # Positions
#             [np.array([0,1,0]),np.array([0,1,0])],                                         # Velocities
#             isRigid=True)
#             ).combine(dumbbell_mesh(
    # [np.array([.5,0.*np.pi/2.-.1,1.*np.pi/2.]),np.array([.5,2.*np.pi/2.-.1,3.*np.pi/2.])],      # Positions
    # [np.array([0,0,0]),np.array([0,-1,0])],isRigid=True)                                         # Velocities
# ))

sim_system = dumbbell_mesh(
    [np.array([.5,np.pi/2.,1*np.pi/2.]),np.array([.5,np.pi/2.,3.*np.pi/2.])],      # Positions
    [np.array([0,0,0]),np.array([0,0,1])],                                         # Velocities
    isRigid=True).combine(dumbbell_mesh(
    [np.array([.5,np.pi/2.,0.*np.pi/2.]),np.array([.5,np.pi/2.,2.*np.pi/2.])],      # Positions
    [np.array([0,0,0]),np.array([0,0,1])],                                         # Velocities
    isRigid=True)).combine(dumbbell_mesh(
    [np.array([.5,0.*np.pi/2.-.1,3.*np.pi/2.]),np.array([.5,2.*np.pi/2.-.1,3.*np.pi/2.])],      # Positions
    [np.array([0,0,0]),np.array([0,-1,0])],isRigid=True)                                         # Velocities
)

# System Parameters
stiffness = 1.
eq_len    = 1.

radii = []
masses = []
for a in range(len(sim_system.data)):
    radius = .2 
    mass   = 1. 
    radii.append(radius)
    masses.append(mass)


# Construct configuration space of system
configurationSpace = ConfigurationSpace(masses, radii, ambientSpace)


dt = .001
tmaxNum = 10000


before_dat = sim_system.clone()

# //make the simulation
sim = Simulation( ambientSpace, sim_system, configurationSpace, dt , "RigidRadau2")

for c in range(tmaxNum):
    if c%100 == 0:
        print(c)
    sim.step()

after_dat = sim_system.clone()

# //make the visualization of the simulation
# let viz = new RenderSim( sim, radii );
    
## Plot Space Trajectory and Distance Error
fig = plt.figure(figsize=(15,6))

## Plot the Space Trajectory of the Rod Body System
ax1=fig.add_subplot(1,3,1,projection='3d')

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

# Data containers for plotting
hyppartdata = [ [] for _ in range(len(sim_system.data)) ]
partdata = [ [] for _ in range(len(sim_system.data)) ]

# hyppart1plot=[]
# hyppart2plot=[]
# hyppart3plot=[]
# hyppart4plot=[]
# part1plot=[]
# part2plot=[]
# part3plot=[]
# part4plot=[]

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
for a in range(len(sim.data_container)):
    for b in range(len(sim.data_container[a].data)):
        hyppartdata[b].append(list(rot2hyp(sim.data_container[a].data[b].pos)))
        partdata[b].append(list(hyp2poin3d(rot2hyp(sim.data_container[a].data[b].pos))))

for c in range(len(partdata)):
    hyppartdata[c] = np.array(hyppartdata[c])
    partdata[c] = np.array(partdata[c])

# for a in sim.data_container:
#     hyppart1plot.append(list(rot2hyp(a.data[0].pos)))
#     hyppart2plot.append(list(rot2hyp(a.data[1].pos)))
#     hyppart3plot.append(list(rot2hyp(a.data[2].pos)))
#     hyppart4plot.append(list(rot2hyp(a.data[3].pos)))

#     part1plot.append(list(hyp2poin3d(rot2hyp(a.data[0].pos))))
#     part2plot.append(list(hyp2poin3d(rot2hyp(a.data[1].pos))))
#     part3plot.append(list(hyp2poin3d(rot2hyp(a.data[2].pos))))
#     part4plot.append(list(hyp2poin3d(rot2hyp(a.data[3].pos))))

# hyppart1plot = np.array(hyppart1plot)
# hyppart2plot = np.array(hyppart2plot)
# hyppart3plot = np.array(hyppart3plot)
# hyppart4plot = np.array(hyppart4plot)
# part1plot = np.array(part1plot)
# part2plot = np.array(part2plot)
# part3plot = np.array(part3plot)
# part4plot = np.array(part4plot)

# Generate the pseudo-shell surface of the vertices of the rod body
shelldata = [ [] for _ in range(len(partdata)) ]

for d in range(len(partdata)):
    partx,party,partz = genballh3(hyppartdata[d][-1],radii[d])
    shelldata[d].append(partx)
    shelldata[d].append(party)
    shelldata[d].append(partz)

# part1x,part1y,part1z=genballh3(hyppart1plot[-1],radii[0])
# part2x,part2y,part2z=genballh3(hyppart2plot[-1],radii[1])
# part3x,part3y,part3z=genballh3(hyppart3plot[-1],radii[2])
# part4x,part4y,part4z=genballh3(hyppart4plot[-1],radii[3])


# # Draw the pseudo-shell surface of the vertices of the rod body and draw trajectories
for g in range(len(partdata)):
    if g == 0 or g == 1:
        ax1.plot_surface(shelldata[g][0], shelldata[g][1], shelldata[g][2], color="b")
    elif g == 2 or g == 3:
        ax1.plot_surface(shelldata[g][0], shelldata[g][1], shelldata[g][2], color="r")
    elif g == 4 or g == 5:
        ax1.plot_surface(shelldata[g][0], shelldata[g][1], shelldata[g][2], color="g")
    ax1.plot3D(partdata[g][:,0], partdata[g][:,1], partdata[g][:,2], label="particle {}".format(g+1))

# ax1.plot_surface(part1x, part1y, part1z, color="b")
# ax1.plot_surface(part2x, part2y, part2z, color="b")
# ax1.plot_surface(part3x, part3y, part3z, color="r")
# ax1.plot_surface(part4x, part4y, part4z, color="r")


# Draw trajectory of vertices of rod body
# ax1.plot3D(part1plot[:,0],part1plot[:,1],part1plot[:,2], label="particle 1")
# ax1.plot3D(part2plot[:,0],part2plot[:,1],part2plot[:,2], label="particle 2")
# ax1.plot3D(part3plot[:,0],part3plot[:,1],part3plot[:,2], label="particle 3")
# ax1.plot3D(part4plot[:,0],part4plot[:,1],part4plot[:,2], label="particle 4")
ax1.legend(loc= 'lower left')

## Plot the error in distance between vertices compared to fixed distance
ax2=fig.add_subplot(1,3,2)
ax3=fig.add_subplot(1,3,3)

# Generate distance error data
condata = [ [] for _ in range(len(sim.data_container[0].rig_connectivity)) ]
dtcondata = [ [] for _ in range(len(sim.data_container[0].rig_connectivity)) ]

for a in range(len(sim.data_container)):
    for b in sim.data_container[a].rig_connectivity:
            # Distance Function
            d12 = ambientSpace.geometry.funcDict["d12"](sim.data_container[a].data[b[0]],sim.data_container[a].data[b[1]])
            # First derivatives of distance function
            da1d12 = ambientSpace.geometry.funcDict["da1d12"](sim.data_container[a].data[b[0]],sim.data_container[a].data[b[1]])
            db1d12 = ambientSpace.geometry.funcDict["db1d12"](sim.data_container[a].data[b[0]],sim.data_container[a].data[b[1]])
            dg1d12 = ambientSpace.geometry.funcDict["dg1d12"](sim.data_container[a].data[b[0]],sim.data_container[a].data[b[1]])
            da2d12 = ambientSpace.geometry.funcDict["da1d12"](sim.data_container[a].data[b[1]],sim.data_container[a].data[b[0]])
            db2d12 = ambientSpace.geometry.funcDict["db1d12"](sim.data_container[a].data[b[1]],sim.data_container[a].data[b[0]])
            dg2d12 = ambientSpace.geometry.funcDict["dg1d12"](sim.data_container[a].data[b[1]],sim.data_container[a].data[b[0]])

            dterms = [da1d12,db1d12,dg1d12,da2d12,db2d12,dg2d12]

            dtd12  = ambientSpace.geometry.funcDict["dtd12"](sim.data_container[a].data[b[0]],sim.data_container[a].data[b[1]],dterms)

            rig_con   = ambientSpace.geometry.funcDict["rig_con"](1, d12)
            dtrig_con = ambientSpace.geometry.funcDict["rig_con_tderivative1"](1, d12, dtd12)

            # print("rig_con")
            # print(rig_con)
            condata[sim.data_container[a].rig_connectivity.index(b)].append(rig_con)
            dtcondata[sim.data_container[a].rig_connectivity.index(b)].append(dtrig_con)

# for a in range(len(sim.data_container)):
#     for b in range(len(sim.data_container[a].rig_connectivity)):
#         # print(b)
#         distdata[b].append(
#             ambientSpace.distance(sim.data_container[a].data[sim.data_container[a].rig_connectivity[b][0]].pos,sim.data_container[a].data[sim.data_container[a].rig_connectivity[b][1]].pos)
#             )
        
# print(distdata)

# counter = 0
# for c in range(len(sim.data_container[a].data)):
#     distdata[counter] = ambientSpace.distance(rot2hyp(data1[b][0:3]),rot2hyp(data1[b][3:6]))
#     counter += 1

# Draw distance error as a function of simulation time
ax2.plot(np.arange(0,tmaxNum*dt+dt,dt),condata[:][0],color='b',label = "Pair 1")
ax2.plot(np.arange(0,tmaxNum*dt+dt,dt),condata[:][1],color='r',label = "Pair 2")
ax2.plot(np.arange(0,tmaxNum*dt+dt,dt),condata[:][2],color='g',label = "Pair 3")
ax2.legend()
ax2.set_title('Constraint Data')
ax2.set_xlabel('t')
ax2.set_ylabel('Position Constraint')

# Draw distance error as a function of simulation time
ax3.plot(np.arange(0,tmaxNum*dt+dt,dt),dtcondata[:][0],color='b',label = "Pair 1")
ax3.plot(np.arange(0,tmaxNum*dt+dt,dt),dtcondata[:][1],color='r',label = "Pair 2")
ax3.plot(np.arange(0,tmaxNum*dt+dt,dt),dtcondata[:][2],color='g',label = "Pair 3")
ax3.legend()
ax3.set_title('Constraint Data')
ax3.set_xlabel('t')
ax3.set_ylabel('Velocity Constraint')

# # Draw distance error as a function of simulation time
# ax2.plot(np.arange(0,tmaxNum*dt+dt,dt),np.full(np.array(distdata[0][:]).shape,np.array(distdata[0][0])) - np.array(distdata[0][:]),color='b',label = "Pair 1")
# ax2.plot(np.arange(0,tmaxNum*dt+dt,dt),np.full(np.array(distdata[1][:]).shape,np.array(distdata[1][0])) - np.array(distdata[1][:]),color='r',label = "Pair 2")
# ax2.plot(np.arange(0,tmaxNum*dt+dt,dt),np.full(np.array(distdata[2][:]).shape,np.array(distdata[2][0])) - np.array(distdata[2][:]),color='g',label = "Pair 3")
# ax2.legend()
# ax2.set_title('Distance Data')
# ax2.set_xlabel('t')
# ax2.set_ylabel('Vertex Separation')

fig.tight_layout()	

plt.show()

# # ------------------------------------------------------------------
# ### Uncomment to just generate gif of trajectory of the particle ###
# # ------------------------------------------------------------------

# # # Generate gif

# # First set up the figure, the axis, and the plot element we want to animate
# fig = plt.figure(figsize=(8,8))
# ax1 = fig.add_subplot(111, projection='3d')
# # ax1.set_aspect("equal")

# # #draw sphere
# # Draw Unit Sphere Horizon in Embedding Space (using Poincare Disk Model of Hyperbolic Space)
# u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
# x = 1.*np.sin(u)*np.cos(v)
# y = 1.*np.sin(u)*np.sin(v)
# z = 1.*np.cos(u)
# ax1.plot_wireframe(x, y, z, color="b", alpha=.1)
# ax1.set_xlim3d(-1,1)
# ax1.set_xlabel('X')
# ax1.set_ylim3d(-1,1)
# ax1.set_ylabel('Y')
# ax1.set_zlim3d(-1,1)
# ax1.set_zlabel('Z')

# # Bounding box
# ub, vb = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
# xb = np.sinh(2)/(np.cosh(2)+1) * np.sin(ub)*np.cos(vb)
# yb = np.sinh(2)/(np.cosh(2)+1) * np.sin(ub)*np.sin(vb)
# zb = np.sinh(2)/(np.cosh(2)+1) * np.cos(ub)
# ax1.plot_wireframe(xb, yb, zb, color="k", alpha=.1)

# # Particle Plot data
# ball_box=[]
# viz_balls=[]

# for h in range(len(partdata)):
#     ball_box.append(genballh3(hyppartdata[h][0],radii[h]))
#     if h == 0 or h == 1:
#         viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")])
#     elif h == 2 or h == 3:
#         viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="r")])
#     elif h == 4 or h == 5:
#         viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="g")])

# # ball_box.append(genballh3(hyppart1plot[0],radii[0])) #hypercirch3(array([sinh(gat[a])*sin(gbt[a])*cos(ggt[a]),sinh(gat[a])*sin(gbt[a])*sin(ggt[a]),sinh(gat[a])*cos(gbt[a]),cosh(gat[a])]),.1))
# # viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")])

# # ball_box.append(genballh3(hyppart2plot[0],radii[1])) #hypercirch3(array([sinh(gat[a])*sin(gbt[a])*cos(ggt[a]),sinh(gat[a])*sin(gbt[a])*sin(ggt[a]),sinh(gat[a])*cos(gbt[a]),cosh(gat[a])]),.1))
# # viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")])

# # ball_box.append(genballh3(hyppart3plot[0],radii[2])) #hypercirch3(array([sinh(gat[a])*sin(gbt[a])*cos(ggt[a]),sinh(gat[a])*sin(gbt[a])*sin(ggt[a]),sinh(gat[a])*cos(gbt[a]),cosh(gat[a])]),.1))
# # viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="r")])

# # ball_box.append(genballh3(hyppart4plot[0],radii[3])) #hypercirch3(array([sinh(gat[a])*sin(gbt[a])*cos(ggt[a]),sinh(gat[a])*sin(gbt[a])*sin(ggt[a]),sinh(gat[a])*cos(gbt[a]),cosh(gat[a])]),.1))
# # viz_balls.append([ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="r")])

# # #draw trajectory
# # for b in range(vert_num):
# #     ax1.plot3D(gut[b::vert_num],gvt[b::vert_num],grt[b::vert_num], label="particle []".format(b))
# # ax1.legend(loc= 'lower left')

# # part1x,part1y,part1z=hypercirch3(array([sinh(gat[0::3][0])*sin(gbt[0::3][0])*cos(ggt[0::3][0]),sinh(gat[0::3][0])*sin(gbt[0::3][0])*sin(ggt[0::3][0]),sinh(gat[0::3][0])*cos(gbt[0::3][0]),cosh(gat[0::3][0])]),particles[0][7])
# # part2x,part2y,part2z=hypercirch3(array([sinh(gat[1::3][1])*sin(gbt[1::3][1])*cos(ggt[1::3][1]),sinh(gat[1::3][1])*sin(gbt[1::3][1])*sin(ggt[1::3][1]),sinh(gat[1::3][1])*cos(gbt[1::3][1]),cosh(gat[1::3][1])]),particles[1][7])
# # part3x,part3y,part3z=hypercirch3(array([sinh(gat[2::3][2])*sin(gbt[2::3][2])*cos(ggt[2::3][2]),sinh(gat[2::3][2])*sin(gbt[2::3][2])*sin(ggt[2::3][2]),sinh(gat[2::3][2])*cos(gbt[2::3][2]),cosh(gat[2::3][2])]),particles[2][7])
# # ball1=[ax1.plot_surface(part1x,part1y,part1z, color="b")]
# # ball2=[ax1.plot_surface(part2x,part2y,part2z, color="r")]
# # ball3=[ax1.plot_surface(part3x,part3y,part3z, color="k")]

# # animation function. This is called sequentially
# frames=100
# def animate(i):
#     for j in range(len(partdata)):
#         ax1.plot3D(partdata[j][:int(tmaxNum*i/frames),0],partdata[j][:int(tmaxNum*i/frames),1],partdata[j][:int(tmaxNum*i/frames),2])

#         ball_box.append(genballh3(hyppartdata[j][int(tmaxNum*i/frames)],radii[j])) #array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
#         viz_balls[j][0].remove()
#         if j == 0 or j == 1:
#             viz_balls[j][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")
#         elif j == 2 or j == 3:
#             viz_balls[j][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="r")
#         elif j == 4 or j == 5:
#             viz_balls[j][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="g")

#     # ax1.plot3D(part1plot[:int(tmaxNum*i/frames),0],part1plot[:int(tmaxNum*i/frames),1],part1plot[:int(tmaxNum*i/frames),2])
#     # ax1.plot3D(part2plot[:int(tmaxNum*i/frames),0],part2plot[:int(tmaxNum*i/frames),1],part2plot[:int(tmaxNum*i/frames),2])
#     # ax1.plot3D(part3plot[:int(tmaxNum*i/frames),0],part3plot[:int(tmaxNum*i/frames),1],part3plot[:int(tmaxNum*i/frames),2])
#     # ax1.plot3D(part4plot[:int(tmaxNum*i/frames),0],part4plot[:int(tmaxNum*i/frames),1],part4plot[:int(tmaxNum*i/frames),2])
#     # ax1.legend(loc= 'lower left')

#     # ball_box.append(genballh3(hyppart1plot[int(tmaxNum*i/frames)],radii[0])) #array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
#     # viz_balls[0][0].remove()
#     # viz_balls[0][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")
#     # ball_box.append(genballh3(hyppart2plot[int(tmaxNum*i/frames)],radii[1])) #array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
#     # viz_balls[1][0].remove()
#     # viz_balls[1][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")
#     # ball_box.append(genballh3(hyppart3plot[int(tmaxNum*i/frames)],radii[0])) #array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
#     # viz_balls[2][0].remove()
#     # viz_balls[2][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="r")
#     # ball_box.append(genballh3(hyppart4plot[int(tmaxNum*i/frames)],radii[1])) #array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
#     # viz_balls[3][0].remove()
#     # viz_balls[3][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="r")
    
#     # ax1.plot3D(gut[b::vert_num][:int(len(timearr)*i/frames)],gvt[b::vert_num][:int(len(timearr)*i/frames)],grt[b::vert_num][:int(len(timearr)*i/frames)], label="particle 1")
#     # for b in range(vert_num):
#     #     ax1.plot3D(gut[b::vert_num][:int(len(timearr)*i/frames)],gvt[b::vert_num][:int(len(timearr)*i/frames)],grt[b::vert_num][:int(len(timearr)*i/frames)], label="particle {}".format(b))
#     # ax1.plot3D(gut[0::3][:int(len(timearr)*i/frames)],gvt[0::3][:int(len(timearr)*i/frames)],grt[0::3][:int(len(timearr)*i/frames)])
#     # ax1.plot3D(gut[1::3][:int(len(timearr)*i/frames)],gvt[1::3][:int(len(timearr)*i/frames)],grt[1::3][:int(len(timearr)*i/frames)])
#     # ax1.plot3D(gut[2::3][:int(len(timearr)*i/frames)],gvt[2::3][:int(len(timearr)*i/frames)],grt[2::3][:int(len(timearr)*i/frames)])
#         # part_plot.append(genballh3(array([sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*cos(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*sin(gbt[b::vert_num][int(len(timearr)*i/frames)])*sin(ggt[b::vert_num][int(len(timearr)*i/frames)]),sinh(gat[b::vert_num][int(len(timearr)*i/frames)])*cos(gbt[b::vert_num][int(len(timearr)*i/frames)]),cosh(gat[b::vert_num][int(len(timearr)*i/frames)])]),.1))
#         # viz_balls[b][0].remove()
#         # viz_balls[b][0]=ax1.plot_surface(ball_box[-1][0], ball_box[-1][1], ball_box[-1][2], color="b")
#     # part1x,part1y,part1z=hypercirch3(array([sinh(gat[0::3][int(len(timearr)*i/frames)])*sin(gbt[0::3][int(len(timearr)*i/frames)])*cos(ggt[0::3][int(len(timearr)*i/frames)]),sinh(gat[0::3][int(len(timearr)*i/frames)])*sin(gbt[0::3][int(len(timearr)*i/frames)])*sin(ggt[0::3][int(len(timearr)*i/frames)]),sinh(gat[0::3][int(len(timearr)*i/frames)])*cos(gbt[0::3][int(len(timearr)*i/frames)]),cosh(gat[0::3][int(len(timearr)*i/frames)])]),particles[0][7])
#     # part2x,part2y,part2z=hypercirch3(array([sinh(gat[1::3][int(len(timearr)*i/frames)])*sin(gbt[1::3][int(len(timearr)*i/frames)])*cos(ggt[1::3][int(len(timearr)*i/frames)]),sinh(gat[1::3][int(len(timearr)*i/frames)])*sin(gbt[1::3][int(len(timearr)*i/frames)])*sin(ggt[1::3][int(len(timearr)*i/frames)]),sinh(gat[1::3][int(len(timearr)*i/frames)])*cos(gbt[1::3][int(len(timearr)*i/frames)]),cosh(gat[1::3][int(len(timearr)*i/frames)])]),particles[1][7])
#     # part3x,part3y,part3z=hypercirch3(array([sinh(gat[2::3][int(len(timearr)*i/frames)])*sin(gbt[2::3][int(len(timearr)*i/frames)])*cos(ggt[2::3][int(len(timearr)*i/frames)]),sinh(gat[2::3][int(len(timearr)*i/frames)])*sin(gbt[2::3][int(len(timearr)*i/frames)])*sin(ggt[2::3][int(len(timearr)*i/frames)]),sinh(gat[2::3][int(len(timearr)*i/frames)])*cos(gbt[2::3][int(len(timearr)*i/frames)]),cosh(gat[2::3][int(len(timearr)*i/frames)])]),particles[2][7])
#     # ball1[0].remove()
#     # ball1[0]=ax1.plot_surface(part1x,part1y,part1z, color="b")
#     # ball2[0].remove()
#     # ball2[0]=ax1.plot_surface(part2x,part2y,part2z, color="r")
#     # ball3[0].remove()
#     # ball3[0]=ax1.plot_surface(part3x,part3y,part3z, color="k")

# # equivalent to rcParams['animation.html'] = 'html5'
# rc('animation', html='html5')

# # call the animator. blit=True means only re-draw the parts that 
# # have changed.
# anim = animation.FuncAnimation(fig, animate,frames=frames, interval=50)

# anim.save('./h3rigidmeshrs2newtest3.mov', writer='imagemagick')