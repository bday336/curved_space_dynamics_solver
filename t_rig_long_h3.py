# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from rod_derivative_bank import *
from triangle_derivative_bank_h3 import *
from triangle_derivative_bank_s3 import *
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3, h3rads2tridae, h3rads3tridae, s3rads2tridae, s3rads3tridae

# Solver Setup

# Time array based on time step
dt = .001    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Time array based on number of steps
# nump = 10000    # Number of steps
# t_max = 10      # Total simulation time
# t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

# Simulation data container

# Sim H3
rs2simdatalist = np.zeros((t_arr.shape[0],18+3))
rs3simdatalist = np.zeros((t_arr.shape[0],18+3))

# Sim S3
# rs2simdatalist = np.zeros((t_arr.shape[0],18+3))
# rs3simdatalist = np.zeros((t_arr.shape[0],18+3))

# Initial Data
v = 1.      # Initial Velocity
# x = 1.      # Spring Rest Length H3
x1 = 1.      # Rod Length H3
x2 = 1.      # Rod Length H3
x3 = 1.      # Rod Length H3
m1 = 1.      # Mass of point masses
m2 = 1.      # Mass of point masses
m3 = 1.      # Mass of point masses
params = [m1,m2,m3,x1,x2,x3]


# Sim bar in H3
startvec = np.array([
    [np.arccosh(np.cosh(x1)/np.cosh(x1/2.)),np.pi/2.,0.],[.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
    killingvech3([np.arccosh(np.cosh(x1)/np.cosh(x1/2.)),np.pi/2.,0.],v,"x"), killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()
startvec = np.append(startvec,0.) #lam1
startvec = np.append(startvec,0.) #lam2
startvec = np.append(startvec,0.) #lam3

# Sim bar in S3
# startvec = np.array([
#     [np.pi/2., np.pi/2., np.arccos(cos(x1)/cos(x1/2.))],[(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
#     killingvecs3([np.pi/2., np.pi/2., np.arccos(cos(x1)/cos(x1/2.))],-v,"vz"), killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()
# startvec = np.append(startvec,0.) #lam1
# startvec = np.append(startvec,0.) #lam2
# startvec = np.append(startvec,0.) #lam3

# Sim in H3
rs2simdatalist[0] = startvec.copy()
rs3simdatalist[0] = startvec.copy()

# Sim in S3
# rs2simdatalist[0] = startvec.copy()
# rs3simdatalist[0] = startvec.copy()


# First Step
step = 1

# Sim in H3
rs2simdatalist[step] = h3rads2tridae(startvec=startvec,params=params,dt=dt)
rs3simdatalist[step] = h3rads3tridae(startvec=startvec,params=params,dt=dt)

# Sim in S3
# rs2simdatalist[step] = s3rads2tridae(startvec=startvec,params=params,dt=dt)
# rs3simdatalist[step] = s3rads3tridae(startvec=startvec,params=params,dt=dt)

# Sim in H3
startvecrs2sim = rs2simdatalist[step]
startvecrs3sim = rs3simdatalist[step]

# Sim in S3
# startvecrs2sim = rs2simdatalist[step]
# startvecrs3sim = rs3simdatalist[step]

step += 1

while (step <= int(t_max/dt)):

    # Sim in H3
    rs2simdatalist[step] = h3rads2tridae(startvec=startvecrs2sim,params=params,dt=dt)
    rs3simdatalist[step] = h3rads3tridae(startvec=startvecrs3sim,params=params,dt=dt)

    # Sim in S3
    # rs2simdatalist[step] = s3rads2tridae(startvec=startvecrs2sim,params=params,dt=dt)
    # rs3simdatalist[step] = s3rads3tridae(startvec=startvecrs3sim,params=params,dt=dt)

    # Sim in H3
    startvecrs2sim = rs2simdatalist[step]
    startvecrs3sim = rs3simdatalist[step]

    # Sim in S3
    # startvecrs2sim = rs2simdatalist[step]
    # startvecrs3sim = rs3simdatalist[step]

    if step%int(1/dt)==0:
            print(step)
    step += 1


# Sim in H3
np.save("h3_rig_t_rads2_sim_tmax10_dt001",rs2simdatalist)
np.save("h3_rig_t_rads3_sim_tmax10_dt001",rs3simdatalist)

# Sim in S3
# np.save("s3_rig_t_rads2_sim_tmax10_dt001",rs2simdatalist)
# np.save("s3_rig_t_rads3_sim_tmax10_dt001",rs3simdatalist)

# data1 = np.load("h3_r_gausss1_sim_tmax10_dt00001.npy")
# data2 = np.load("h3_r_gausss2_sim_tmax10_dt00001.npy")
# data3 = np.load("h3_r_gausss3_sim_tmax10_dt00001.npy")

# data1 = np.load("h3_r_gausss3_gt_tmax10_dt000005.npy")

data2 = np.load("h3_rig_t_rads2_sim_tmax10_dt001.npy")
data3 = np.load("h3_rig_t_rads3_sim_tmax10_dt001.npy")

# data2 = np.load("s3_rig_t_rads2_sim_tmax10_dt001.npy")
# data3 = np.load("s3_rig_t_rads3_sim_tmax10_dt001.npy")

fig,ax=plt.subplots(1,1)


distdatasp12r2 = np.zeros(np.shape(data2)[0])
distdatasp12r3 = np.zeros(np.shape(data3)[0])

distdatasp13r2 = np.zeros(np.shape(data2)[0])
distdatasp13r3 = np.zeros(np.shape(data3)[0])

distdatasp23r2 = np.zeros(np.shape(data2)[0])
distdatasp23r3 = np.zeros(np.shape(data3)[0])


# counter = 0
# for a in range(np.shape(data2)[0]):
#     distdatasp12r2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))-1.
#     distdatasp12r3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))-1.

#     distdatasp13r2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][6:9]))-1.
#     distdatasp13r3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][6:9]))-1.

#     distdatasp23r2[counter] = r4dist(rot2r4(data2[a][3:6]),rot2r4(data2[a][6:9]))-1.
#     distdatasp23r3[counter] = r4dist(rot2r4(data3[a][3:6]),rot2r4(data3[a][6:9]))-1.
#     counter += 1

counter = 0
for a in range(np.shape(data2)[0]):
    distdatasp12r2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))-1.
    distdatasp12r3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))-1.

    distdatasp13r2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][6:9]))-1.
    distdatasp13r3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][6:9]))-1.

    distdatasp23r2[counter] = h3dist(rot2hyp(data2[a][3:6]),rot2hyp(data2[a][6:9]))-1.
    distdatasp23r3[counter] = h3dist(rot2hyp(data3[a][3:6]),rot2hyp(data3[a][6:9]))-1.
    counter += 1


# ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")
# ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,distdata1,'r',label = "Gauss s1")

ax.plot(t_arr,distdatasp12r2,marker = ".",color='k',linestyle = "None",label = "Radau s2 d12")
ax.plot(t_arr,distdatasp12r3,marker = ".",color='b',linestyle = "None",label = "Radau s3 d12")
ax.plot(t_arr,distdatasp13r2,marker = ".",color='k',linestyle = "None",label = "Radau s2 d13")
ax.plot(t_arr,distdatasp13r3,marker = ".",color='b',linestyle = "None",label = "Radau s3 d13")
ax.plot(t_arr,distdatasp23r2,marker = ".",color='k',linestyle = "None",label = "Radau s2 d23")
ax.plot(t_arr,distdatasp23r3,marker = ".",color='b',linestyle = "None",label = "Radau s3 d23")

ax.legend()
ax.set_title('Simulation Data')
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







