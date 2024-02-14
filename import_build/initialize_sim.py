
# import { randomVec3Ball } from "./utils/random.js";
# import {State} from "./Computation/State.js";

# import {DataList} from "./Computation/DataList.js";

# import {ConfigurationSpace} from "./ConfigurationSpace/ConfigurationSpace.js";
# import {Simulation} from "./ConfigurationSpace/Simulation.js";
# import {RenderSim} from "./Visualization/RenderSim.js";

# import {euclidean} from "./AmbientSpace/ExampleSpaces/Euclidean.js";
# import { hyperbolic } from "./AmbientSpace/ExampleSpaces/HypSpacePoincareBall.js";
# import { spherical } from "./AmbientSpace/ExampleSpaces/SphericalStereoProj.js";
# import { inhomogeneousNeg } from "./AmbientSpace/ExampleSpaces/InhomogeneousNeg.js";
# import { inhomogeneousPos } from "./AmbientSpace/ExampleSpaces/InhomogeneousPos.js";
# import { h2xe } from "./AmbientSpace/ExampleSpaces/H2xE.js";
# import { s2xe } from "./AmbientSpace/ExampleSpaces/S2xE-Stereo.js";

# Allow for package importing
import os, sys
sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.path.dirname(os.getcwd())+"/src")

def genball(center,rad):
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

import numpy as np
import matplotlib.pyplot as plt
from random import random,uniform

# Import Packages
from src.Computation.DataList import DataList

from src.Computation.State import State

from src.ConfigurationSpace.ConfigurationSpace import ConfigurationSpace
from src.ConfigurationSpace.ConfigurationSpace import ConfigurationSpace
from src.ConfigurationSpace.Simulation import Simulation

from src.AmbientSpace.ExampleSpaces.Euclidean import euclidean

# Set the ambient space for the simulation environment
ambientSpace = euclidean

# Build a configuration space
NumBalls = 2
MaxRad = 2 #ambientSpace.obstacle.size/5.

radii = []
masses = []
for a in range(NumBalls):
    r = MaxRad #* random()+0.05
    m = 1. #10.*r*r*r
    radii.append(r)
    masses.append(m)


configurationSpace = ConfigurationSpace(masses, radii, ambientSpace)
# // let maxPos = 2.;
# // let maxVel = 1;

# //build the initial set of states for the system:
iniCond = [
    State(np.array([4,0,0]),np.array([-1,0,0])),
    State(np.array([0,4,1]),np.array([0,-1,0]))
]

dataList = DataList(iniCond)

# //make the simulation
sim = Simulation( ambientSpace, configurationSpace, dataList, 0.001 )

for c in range(10000):
    sim.step()

# //make the visualization of the simulation
# let viz = new RenderSim( sim, radii );
    
## Plot Space Trajectory and Distance Error
fig = plt.figure(figsize=(6,6))

## Plot the Space Trajectory of the Rod Body System
ax1=fig.add_subplot(1,1,1,projection='3d')

# Draw Unit Sphere Horizon in Embedding Space (using Poincare Disk Model of Hyperbolic Space)
u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
x = 6.*np.sin(u)*np.cos(v)
y = 6.*np.sin(u)*np.sin(v)
z = 6.*np.cos(u)
ax1.plot_wireframe(x, y, z, color="b", alpha=.1)
ax1.set_xlim3d(-10,10)
ax1.set_xlabel('X')
ax1.set_ylim3d(-10,10)
ax1.set_ylabel('Y')
ax1.set_zlim3d(-10,10)
ax1.set_zlabel('Z')

part1plot=[]
part2plot=[]

for a in sim.data_container:
    part1plot.append(a.data[0].pos)
    part2plot.append(a.data[1].pos)

part1plot = np.array(part1plot)
part2plot = np.array(part2plot)

# Generate the pseudo-shell surface of the vertices of the rod body
part1x,part1y,part1z=genball(part1plot[-1,0:3],radii[0])
part2x,part2y,part2z=genball(part2plot[-1,0:3],radii[1])


# Draw the pseudo-shell surface of the vertices of the rod body
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

# //send the visualization off to be rendered on the screen
# export default { viz };


# //export these to use as global variables in the DEFINITIONS of the classes
# //Simulation, ConfigurationSpace, and RenderSim :0 :0 PLZ FIX
# export { ambientSpace, configurationSpace };