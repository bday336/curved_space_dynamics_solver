# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, hypercirch3, rot2r4, r42rot, s2rstproj, r22s2stproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3simcenpot, dynjac_h3simcenpot, dynfunc_s3simcenpot, dynjac_s3simcenpot

# Solver Setup

# Time array based on time step
dt = .01    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Time array based on number of steps
# nump = 10000    # Number of steps
# t_max = 10      # Total simulation time
# t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

# Simulation data container

# Sim H3
# gs1simdatalist = np.zeros((t_arr.shape[0],6))
# gs2simdatalist = np.zeros((t_arr.shape[0],6))
# gs3simdatalist = np.zeros((t_arr.shape[0],6))

# Sim S3
gs1simdatalist = np.zeros((t_arr.shape[0],6))
gs2simdatalist = np.zeros((t_arr.shape[0],6))
gs3simdatalist = np.zeros((t_arr.shape[0],6))

# Initial Data
v = 1.2     # Initial Velocity
G = 1.     # Newton Gravity constant
ms = 10.    # Mass of potential source
m = 1.     # Mass of test mass
params = [v,G,ms,m]

# Sim bar in H3
# startvec = np.array([
#     [1.,np.pi/2.,np.pi/2.],
#     killingvech3([1.,np.pi/2.,np.pi/2.],v,"x")]).flatten()
# Sim bar in S3 (v=3.5 is closed orbit - they all are as expected >.> from the kepler paper)
# startvec = np.array([
#     [(np.pi-1.)/2.,np.pi/2.,0.],
#     killingvecs3([np.pi/2.,np.pi/2.,0.],-v,"vz")]).flatten()
# S3 test for fig8
startvec = np.array([
    [(np.pi)/2.,np.pi/2.,0.],
    [-2.75,0,3.5]]).flatten()


# Sim in H3
# gs1simdatalist[0] = startvec.copy()
# gs2simdatalist[0] = startvec.copy()
# gs3simdatalist[0] = startvec.copy()

# Sim in S3
gs1simdatalist[0] = startvec.copy()
gs2simdatalist[0] = startvec.copy()
gs3simdatalist[0] = startvec.copy()

# First Step
step = 1

# Sim in H3
# gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3simcenpot,dynjac=dynjac_h3simcenpot,dt=dt)
# gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3simcenpot,dynjac=dynjac_h3simcenpot,dt=dt)
# gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3simcenpot,dynjac=dynjac_h3simcenpot,dt=dt) 

# Sim in S3
gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3simcenpot,dynjac=dynjac_s3simcenpot,dt=dt)
gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3simcenpot,dynjac=dynjac_s3simcenpot,dt=dt)
gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3simcenpot,dynjac=dynjac_s3simcenpot,dt=dt) 

# Sim in H3
# startvecgs1sim = gs1simdatalist[step]
# startvecgs2sim = gs2simdatalist[step]
# startvecgs3sim = gs3simdatalist[step]

# Sim in S3
startvecgs1sim = gs1simdatalist[step]
startvecgs2sim = gs2simdatalist[step]
startvecgs3sim = gs3simdatalist[step]

step += 1

while (step <= t_max/dt):
    # Sim in H3
    # gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_h3simcenpot,dynjac=dynjac_h3simcenpot,dt=dt)
    # gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_h3simcenpot,dynjac=dynjac_h3simcenpot,dt=dt)
    # gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_h3simcenpot,dynjac=dynjac_h3simcenpot,dt=dt)    

    # Sim in S3
    gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_s3simcenpot,dynjac=dynjac_s3simcenpot,dt=dt)
    gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_s3simcenpot,dynjac=dynjac_s3simcenpot,dt=dt)
    gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_s3simcenpot,dynjac=dynjac_s3simcenpot,dt=dt) 

    # Sim in H3
    # startvecgs1sim = gs1simdatalist[step]
    # startvecgs2sim = gs2simdatalist[step]
    # startvecgs3sim = gs3simdatalist[step]

    # Sim in S3
    startvecgs1sim = gs1simdatalist[step]
    startvecgs2sim = gs2simdatalist[step]
    startvecgs3sim = gs3simdatalist[step]

    if step%100==0:
            print(step)
    step += 1

# Sim in H3
# np.save("h3_cenpot_gausss1_sim_tmax10_dt01",gs1simdatalist)
# np.save("h3_cenpot_gausss2_sim_tmax10_dt01",gs2simdatalist)
# np.save("h3_cenpot_gausss3_sim_tmax10_dt01",gs3simdatalist)

# Sim in S3
np.save("s3_cenpot_gausss1_sim_tmax10_dt01",gs1simdatalist)
np.save("s3_cenpot_gausss2_sim_tmax10_dt01",gs2simdatalist)
np.save("s3_cenpot_gausss3_sim_tmax10_dt01",gs3simdatalist)

# data1 = np.load("h3_cenpot_gausss1_sim_tmax10_dt01.npy")
# data2 = np.load("h3_cenpot_gausss2_sim_tmax10_dt01.npy")
# data3 = np.load("h3_cenpot_gausss3_sim_tmax10_dt01.npy")

data1 = np.load("s3_cenpot_gausss1_sim_tmax10_dt01.npy")
data2 = np.load("s3_cenpot_gausss2_sim_tmax10_dt01.npy")
data3 = np.load("s3_cenpot_gausss3_sim_tmax10_dt01.npy")

visdata1 = np.zeros(np.shape(data1[:,0:4]))
visdata2 = np.zeros(np.shape(data2[:,0:4]))
visdata3 = np.zeros(np.shape(data3[:,0:4]))

# H3
# plotdata1 = np.zeros(np.shape(data1[:,0:3]))
# plotdata2 = np.zeros(np.shape(data2[:,0:3]))
# plotdata3 = np.zeros(np.shape(data3[:,0:3]))

# S3
r2plotdata1 = np.zeros(np.shape(data1[:,0:3]))
r2plotdata2 = np.zeros(np.shape(data2[:,0:3]))
r2plotdata3 = np.zeros(np.shape(data3[:,0:3]))
s2plotdata1 = np.zeros(np.shape(data1[:,0:3]))
s2plotdata2 = np.zeros(np.shape(data2[:,0:3]))
s2plotdata3 = np.zeros(np.shape(data3[:,0:3]))

# H3
# counter = 0
# for a in range(np.shape(data1)[0]):
#     visdata1[counter] = rot2hyp(data1[a,0:3])
#     visdata2[counter] = rot2hyp(data2[a,0:3])
#     visdata3[counter] = rot2hyp(data3[a,0:3])
#     plotdata1[counter] = hyp2poin3d(rot2hyp(data1[a,0:3]))
#     plotdata2[counter] = hyp2poin3d(rot2hyp(data2[a,0:3]))
#     plotdata3[counter] = hyp2poin3d(rot2hyp(data3[a,0:3]))
#     counter += 1

# S3
counter = 0
for a in range(np.shape(data1)[0]):
    visdata1[counter] = rot2r4(data1[a,0:3])
    visdata2[counter] = rot2r4(data2[a,0:3])
    visdata3[counter] = rot2r4(data3[a,0:3])
    r2plotdata1[counter] = s2rstproj(rot2r4(data1[a,0:3]))
    r2plotdata2[counter] = s2rstproj(rot2r4(data2[a,0:3]))
    r2plotdata3[counter] = s2rstproj(rot2r4(data3[a,0:3]))
    s2plotdata1[counter] = r22s2stproj(s2rstproj(rot2r4(data1[a,0:3]))[0:2])
    s2plotdata2[counter] = r22s2stproj(s2rstproj(rot2r4(data2[a,0:3]))[0:2])
    s2plotdata3[counter] = r22s2stproj(s2rstproj(rot2r4(data3[a,0:3]))[0:2])
    counter += 1

# H3
# #Plot
# fig = plt.figure(figsize=(8,4))
# ax1 = fig.add_subplot(121, projection='3d')
# # ax1.set_aspect("equal")

# #draw sphere
# u, v = np.mgrid[0:np.pi+(np.pi)/15.:(np.pi)/15., 0:2.*np.pi+(2.*np.pi)/15.:(2.*np.pi)/15.]
# x = np.sin(u)*np.cos(v)
# y = np.sin(u)*np.sin(v)
# z = np.cos(u)
# ax1.plot_wireframe(x, y, z, color="b", alpha=.1)
# ax1.set_xlim3d(-1,1)
# ax1.set_xlabel('X')
# ax1.set_ylim3d(-1,1)
# ax1.set_ylabel('Y')
# ax1.set_zlim3d(-1,1)
# ax1.set_zlabel('Z')

# # If want to plot the test mass sphere
# # part1x,part1y,part1z=hypercirch3(visdata1[-1],.2)
# # part2x,part2y,part2z=hypercirch3(visdata2[-1],.2)
# # part3x,part3y,part3z=hypercirch3(visdata3[-1],.2)

# #draw trajectory
# ax1.plot3D(plotdata1[:,0],plotdata1[:,1],plotdata1[:,2], label="gauss 1", color="b")
# ax1.plot3D(plotdata2[:,0],plotdata2[:,1],plotdata2[:,2], label="gauss 2", color="r")
# ax1.plot3D(plotdata3[:,0],plotdata3[:,1],plotdata3[:,2], label="gauss 3", color="k")
# ax1.legend(loc= 'lower left')

# # If want to plot the test mass sphere
# # ax1.plot_surface(part1x, part1y, part1z, color="b")
# # ax1.plot_surface(part2x, part2y, part2z, color="r")
# # ax1.plot_surface(part3x, part3y, part3z, color="k")

# # Displacement Plot
# ax2=fig.add_subplot(1,2,2)

# ax2.plot(plotdata1[:,0],plotdata1[:,1],'b',label = "Gauss s1")
# ax2.plot(plotdata2[:,0],plotdata2[:,1],'r',label = "Gauss s2")
# ax2.plot(plotdata3[:,0],plotdata3[:,1],'k',label = "Gauss s3")
# ax2.plot(np.cos(np.linspace(0,2*np.pi,100)),np.sin(np.linspace(0,2*np.pi,100)))
# #ax2.axhline(y=spring[3]+sqrt(2.*1./spring[2]*.5*.5), color='b', linestyle='-')
# #ax2.axhline(y=((dist.max()-spring[3])+(dist.min()-spring[3]))/2., color='r', linestyle='-')
# #ax2.set_yscale("log",basey=10)	
# ax2.set_ylabel('y')
# ax2.set_xlabel('x')
# ax2.set_xlim(-1,1)
# ax2.set_ylim(-1,1)
# ax2.legend(loc='lower right')

# plt.show()

# S3
#Plot
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(121, projection='3d')
# ax1.set_aspect("equal")

#draw sphere
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

#draw trajectory
ax1.plot3D(s2plotdata1[:,0],s2plotdata1[:,1],s2plotdata1[:,2], label="gauss 1", color="b")
ax1.plot3D(s2plotdata2[:,0],s2plotdata2[:,1],s2plotdata2[:,2], label="gauss 2", color="r")
ax1.plot3D(s2plotdata3[:,0],s2plotdata3[:,1],s2plotdata3[:,2], label="gauss 3", color="k")
ax1.legend(loc= 'lower left')

# Displacement Plot
ax2=fig.add_subplot(1,2,2)

ax2.plot(s2plotdata1[:,0],s2plotdata1[:,1],'b',label = "Gauss s1")
ax2.plot(s2plotdata2[:,0],s2plotdata2[:,1],'r',label = "Gauss s2")
ax2.plot(s2plotdata3[:,0],s2plotdata3[:,1],'k',label = "Gauss s3")
#ax2.axhline(y=spring[3]+sqrt(2.*1./spring[2]*.5*.5), color='b', linestyle='-')
#ax2.axhline(y=((dist.max()-spring[3])+(dist.min()-spring[3]))/2., color='r', linestyle='-')
#ax2.set_yscale("log",basey=10)	
#ax2.set_ylabel('displacement (m)')
ax2.set_xlabel('time (s)')
ax2.legend(loc='lower right')

plt.show()


# Distance vs time plots
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

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
#     distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
#     distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
#     counter += 1

# ax.plot(t_arr,data1[:,0],'r',label = "Gauss s1")
# ax.plot(t_arr,data2[:,0],'k',label = "Gauss s2")
# ax.plot(t_arr,data3[:,0],'b',label = "Gauss s3")
# # ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
# # ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# # ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()







