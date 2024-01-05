# Allow for package import
import os, sys
sys.path.append(os.path.dirname(os.getcwd()))
sys.path.append(os.path.dirname(os.getcwd())+"/src")

# Import Packages
import numpy as np
import matplotlib.pyplot as plt
from src.function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, hypercirch3
from src.integrator_files.integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from src.test_system_simulations.test_system_bank import dynfunc_h3simbar, dynjac_h3simbar, dynfunc_s3simbar, dynjac_s3simbar
from src.test_system_simulations.SpringBarSimulation import SpringBarSimulation

### Example simulation of a rod body with rigid connection between vertices in 3-dimensional hyperbolic space (H3)
### Using the rotational parameterization of hyperboloid model of H3 embedded in E^(3,1) (4D Minkowski Space)
###     x = sinh(a) * sin(b) * cos(g) 
###     y = sinh(a) * sin(b) * sin(g)
###     z = sinh(a) * cos(b)
###     w = cosh(a)


## Time array based on time step
dt = .1          # Time step size
t_max = 10.      # Total simulation time

## System Data
geometry = "h3"     # Geometry of ambient space (currently set to 3-dimensional hyperbolic space)
v = 1.              # Initial Velocity (perpendicular to connecting geodesic between vertices)
ks = 1.             # Spring Stiffness
x = 1.              # Spring Rest Length
m = 1.              # Mass of point mass (vertex) - same for each
r1 = .2             # Radius of vertex 1
r2 = .2             # Radius of vertex 2
params = [v,ks,x,m]


## Initial Data for System
# Data here is provided in terms of rotational parameterization of H3 (see build.src.function_bank for details)
startvec = np.array([
    [.5,np.pi/2.,np.pi/2.],                         # Initial Position of vertex 1
    [.5,np.pi/2.,3.*np.pi/2.],                      # Initial Position of vertex 2
    killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"),     # Initial Velocity of vertex 1
    killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")   # Initial Velocity of vertex 2
    ]).flatten()


## Solver 
solver_id = "gs1"   # Here using Gauss collocation method with 1 internal step

## Initialize Simulation Object and Run Simulation
sim_test = SpringBarSimulation(geometry, params, dt, t_max, solver_id)
sim_test.set_initial_conditions(startvec)
sim_test.run()
sim_test.output_data()


## Read-In Data File for Data Analysis and Visualization
data1 = np.load("h3_r_{}_sim_tmax{}_dt{}.npy".format(solver_id, str(t_max), str(dt)))

# --------------------------------------------------------------------
### Plot trajectory in the Poincare disk model with distance plots ###
# --------------------------------------------------------------------

## Plot Space Trajectory and Distance Error
fig = plt.figure(figsize=(12,6))

## Plot the Space Trajectory of the Rod Body System
ax1=fig.add_subplot(1,2,1,projection='3d')

# Draw Unit Sphere Horizon in Embedding Space (using Poincare Disk Model of Hyperbolic Space)
u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
x = np.sin(u)*np.cos(v)
y = np.sin(u)*np.sin(v)
z = np.cos(u)
ax1.plot_wireframe(x, y, z, color="b", alpha=.1)
ax1.set_xlim3d(-1,1)
ax1.set_xlabel('X')
ax1.set_ylim3d(-1,1)
ax1.set_ylabel('Y')
ax1.set_zlim3d(-1,1)
ax1.set_zlabel('Z')

# Generate the pseudo-shell surface of the vertices of the rod body
part1x,part1y,part1z=hypercirch3(rot2hyp(data1[-1,:3]),.2)
part2x,part2y,part2z=hypercirch3(rot2hyp(data1[-1,3:6]),.2)

# Transform from Hyperboloid model into Poincare disk model for plotting
part1plot=[]
part2plot=[]

for a in range(data1.shape[0]):
    part1plot.append(list(hyp2poin3d(rot2hyp(data1[a,:3]))))
    part2plot.append(list(hyp2poin3d(rot2hyp(data1[a,3:6]))))

part1plot = np.array(part1plot)
part2plot = np.array(part2plot)

# Draw trajectory of vertices of rod body
ax1.plot3D(part1plot[:,0],part1plot[:,1],part1plot[:,2], label="particle 1")
ax1.plot3D(part2plot[:,0],part2plot[:,1],part2plot[:,2], label="particle 2")
ax1.legend(loc= 'lower left')

# Draw the pseudo-shell surface of the vertices of the rod body
ax1.plot_surface(part1x, part1y, part1z, color="b")
ax1.plot_surface(part2x, part2y, part2z, color="b")

## Plot the error in distance between vertices compared to fixed distance
ax2=fig.add_subplot(1,2,2)

# Generate distance error data
distdata1 = np.zeros(np.shape(data1)[0])

counter = 0
for b in range(np.shape(data1)[0]):
    distdata1[counter] = h3dist(rot2hyp(data1[b][0:3]),rot2hyp(data1[b][3:6]))
    counter += 1

# Draw distance error as a function of simulation time
ax2.plot(sim_test.t_arr,distdata1,marker = ".",color='k',linestyle = "None",label = "Radau s2")
ax2.legend()
ax2.set_title('Simulation Data')
ax2.set_xlabel('t')
ax2.set_ylabel('l')

fig.tight_layout()	

plt.show()

    