# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3simgeo, dynjac_h3simgeo, dynfunc_s3simgeo, dynjac_s3simgeo


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
gs1simdatalist = np.zeros((t_arr.shape[0],6))
gs2simdatalist = np.zeros((t_arr.shape[0],6))
gs3simdatalist = np.zeros((t_arr.shape[0],6))

# Sim S3
# gs1simdatalist = np.zeros((t_arr.shape[0],6))
# gs2simdatalist = np.zeros((t_arr.shape[0],6))
# gs3simdatalist = np.zeros((t_arr.shape[0],6))

# Initial Data
v = 1.      # Initial Velocity
m = 1.      # Mass of point masses
params = [v,m]

# Sim bar in H3
# startvec = np.array([
    # [.5,np.pi/2.,0.],
    # killingvech3([.5,np.pi/2.,0.],v,"x")]).flatten()
# startvec = np.array([
#     [.5,np.pi/2.,np.pi/2.],
#     killingvech3([.5,np.pi/2.,np.pi/2.],v,"x")]).flatten()
# startvec = np.array([
#     [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()

# Sim bar in S3
startvec = np.array([
    [(np.pi-1.)/2.,np.pi/2.,0.],
    killingvecs3([(np.pi-1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()
# startvec = np.array([
#     [(np.pi - 1.)/2.,np.pi/2.,0.],
#     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()
# startvec = np.array([
#     [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
#     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()


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
# gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3simgeo,dynjac=dynjac_h3simgeo,dt=dt,tol=1e-15)
# gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3simgeo,dynjac=dynjac_h3simgeo,dt=dt,tol=1e-15)
# gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3simgeo,dynjac=dynjac_h3simgeo,dt=dt,tol=1e-15)


# Sim in S3
gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3simgeo,dynjac=dynjac_s3simgeo,dt=dt,tol=1e-15)
gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3simgeo,dynjac=dynjac_s3simgeo,dt=dt,tol=1e-15)
gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3simgeo,dynjac=dynjac_s3simgeo,dt=dt,tol=1e-15)


# Sim in H3
# startvecgs1sim = gs1simdatalist[step]
# startvecgs2sim = gs2simdatalist[step]
# startvecgs3sim = gs3simdatalist[step]


# Sim in S3
startvecgs1sim = gs1simdatalist[step]
startvecgs2sim = gs2simdatalist[step]
startvecgs3sim = gs3simdatalist[step]

step += 1

while (step <= int(t_max/dt)):

    # Sim in H3
    # gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_h3simgeo,dynjac=dynjac_h3simgeo,dt=dt,tol=1e-15)
    # gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_h3simgeo,dynjac=dynjac_h3simgeo,dt=dt,tol=1e-15)
    # gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_h3simgeo,dynjac=dynjac_h3simgeo,dt=dt,tol=1e-15)

    # Sim in S3
    gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_s3simgeo,dynjac=dynjac_s3simgeo,dt=dt,tol=1e-15)
    gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_s3simgeo,dynjac=dynjac_s3simgeo,dt=dt,tol=1e-15)
    gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_s3simgeo,dynjac=dynjac_s3simgeo,dt=dt,tol=1e-15)

    # Sim in H3
    # startvecgs1sim = gs1simdatalist[step]
    # startvecgs2sim = gs2simdatalist[step]
    # startvecgs3sim = gs3simdatalist[step]

    # Sim in S3
    startvecgs1sim = gs1simdatalist[step]
    startvecgs2sim = gs2simdatalist[step]
    startvecgs3sim = gs3simdatalist[step]

    if step%int(1/dt)==0:
            print(step)
    step += 1


# Sim in H3
# np.save("h3_geo_gausss1_sim_tmax10_dt01",gs1simdatalist)
# np.save("h3_geo_gausss2_sim_tmax10_dt01",gs2simdatalist)
# np.save("h3_geo_gausss3_sim_tmax10_dt01",gs3simdatalist)

# Sim in S3
np.save("s3_geo_gausss1_sim_tmax10_dt01",gs1simdatalist)
np.save("s3_geo_gausss2_sim_tmax10_dt01",gs2simdatalist)
np.save("s3_geo_gausss3_sim_tmax10_dt01",gs3simdatalist)

# data1 = np.load("h3_geo_gausss1_sim_tmax10_dt00001.npy")
# data2 = np.load("h3_geo_gausss2_sim_tmax10_dt00001.npy")
# data3 = np.load("h3_geo_gausss3_sim_tmax10_dt00001.npy")

# data1 = np.load("s3_r_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("s3_r_gausss2_sim_tmax10_dt1.npy")
# data3 = np.load("s3_r_gausss3_sim_tmax10_dt1.npy")

# data1 = np.load("s3_r_gausss3_gt_tmax10_dt00001.npy")


fig,ax=plt.subplots(1,1)

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

# ax.plot(t_arr,2.*(np.pi/2. - data1[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - data2[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - data3[:,0]),'b',label = "Gauss s3")
# ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,2.*(data2[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,2.*(data3[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
# ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()







