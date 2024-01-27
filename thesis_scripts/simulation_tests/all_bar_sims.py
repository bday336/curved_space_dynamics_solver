# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3exactbar, dynjac_h3exactbar, dynfunc_h3simbar, dynjac_h3simbar, dynfunc_s3exactbar, dynjac_s3exactbar, dynfunc_s3simbar, dynjac_s3simbar

# Solver Setup

# Time array based on time step
dt = .1    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Time array based on number of steps
# nump = 10000    # Number of steps
# t_max = 10      # Total simulation time
# t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

# Simulation data container

# Exact H3
# gs1exactdatalist = np.zeros((t_arr.shape[0],2))
# gs2exactdatalist = np.zeros((t_arr.shape[0],2))
# gs3exactdatalist = np.zeros((t_arr.shape[0],2))

# Sim H3
gs1simdatalist = np.zeros((t_arr.shape[0],12))
gs2simdatalist = np.zeros((t_arr.shape[0],12))
gs3simdatalist = np.zeros((t_arr.shape[0],12))

# Exact S3
# gs1exactdatalist = np.zeros((t_arr.shape[0],2))
# gs2exactdatalist = np.zeros((t_arr.shape[0],2))
# gs3exactdatalist = np.zeros((t_arr.shape[0],2))

# Sim S3
# gs1simdatalist = np.zeros((t_arr.shape[0],12))
# gs2simdatalist = np.zeros((t_arr.shape[0],12))
# gs3simdatalist = np.zeros((t_arr.shape[0],12))

# Initial Data
v = 1.      # Initial Velocity
ks = 1.     # Spring Stiffness
x = 1.      # Spring Rest Length H3
# x = np.pi - 1.      # Spring Rest Length S3 exact
# x = 1.      # Spring Rest Length S3
m = 1.      # Mass of point masses
params = [v,ks,x,m]

# Exact bar in H3
# startvec = np.array([.5,0.])
# Sim bar in H3
# startvec = np.array([
#     [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()

# Exact bar in S3
# startvec = np.array([(np.pi - 1.)/2.,0.])
# Sim bar in S3
startvec = np.array([
    [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
    killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()

# Exact in H3
# gs1exactdatalist[0] = startvec.copy()
# gs2exactdatalist[0] = startvec.copy()
# gs3exactdatalist[0] = startvec.copy()

# Sim in H3
# gs1simdatalist[0] = startvec.copy()
# gs2simdatalist[0] = startvec.copy()
# gs3simdatalist[0] = startvec.copy()

# Exact in S3
# gs1exactdatalist[0] = startvec.copy()
# gs2exactdatalist[0] = startvec.copy()
# gs3exactdatalist[0] = startvec.copy()

# Sim in S3
gs1simdatalist[0] = startvec.copy()
gs2simdatalist[0] = startvec.copy()
gs3simdatalist[0] = startvec.copy()

# First Step
step = 1

# Exact in H3
# gs1exactdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt,tol=1e-15,imax=1e6)
# gs2exactdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt,tol=1e-15,imax=1e6)
# gs3exactdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt,tol=1e-15,imax=1e6)

# Sim in H3
# gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt,tol=1e-15)
# gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt,tol=1e-15)
# gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt,tol=1e-15)

# Exact in S3
# gs1exactdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt,tol=1e-15,imax=1e6)
# gs2exactdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt,tol=1e-15,imax=1e6)
# gs3exactdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt,tol=1e-15,imax=1e6)

# Sim in S3
gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt,tol=1e-15)
gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt,tol=1e-15)
gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt,tol=1e-15)

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
startvecgs1sim = gs1simdatalist[step]
startvecgs2sim = gs2simdatalist[step]
startvecgs3sim = gs3simdatalist[step]

step += 1

while (step <= int(t_max/dt)):
    # Exact in H3
    # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt,tol=1e-15,imax=1e6)
    # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt,tol=1e-15,imax=1e6)
    # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt,tol=1e-15,imax=1e6)

    # Sim in H3
    # gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt,tol=1e-15)
    # gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt,tol=1e-15)
    # gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt,tol=1e-15)

    # Exact in S3
    # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt,tol=1e-15,imax=1e6)
    # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt,tol=1e-15,imax=1e6)
    # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt,tol=1e-15,imax=1e6)

    # Sim in S3
    gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt,tol=1e-15)
    gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt,tol=1e-15)
    gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt,tol=1e-15)

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
    startvecgs1sim = gs1simdatalist[step]
    startvecgs2sim = gs2simdatalist[step]
    startvecgs3sim = gs3simdatalist[step]

    if step%int(1/dt)==0:
            print(step)
    step += 1

# Exact in H3
# np.save("h3_r_gausss1_anal_tmax10_dt00001",gs1exactdatalist)
# np.save("h3_r_gausss2_anal_tmax10_dt00001",gs2exactdatalist)
# np.save("h3_r_gausss3_anal_tmax10_dt00001",gs3exactdatalist)

# Sim in H3
# np.save("h3_r_gausss1_sim_tmax10_dt00001",gs1simdatalist)
# np.save("h3_r_gausss2_sim_tmax10_dt00001",gs2simdatalist)
# np.save("h3_r_gausss3_sim_tmax10_dt00001",gs3simdatalist)

# Exact in S3
# np.save("s3_r_gausss1_anal_tmax10_dt00001",gs1exactdatalist)
# np.save("s3_r_gausss2_anal_tmax10_dt00001",gs2exactdatalist)
# np.save("s3_r_gausss3_anal_tmax10_dt00001",gs3exactdatalist)

# Sim in S3
np.save("s3_r_gausss1_sim_tmax10_dt1",gs1simdatalist)
np.save("s3_r_gausss2_sim_tmax10_dt1",gs2simdatalist)
np.save("s3_r_gausss3_sim_tmax10_dt1",gs3simdatalist)

# data1 = np.load("h3_r_gausss1_anal_tmax10_dt00001.npy")
# data2 = np.load("h3_r_gausss2_anal_tmax10_dt00001.npy")
# data3 = np.load("h3_r_gausss3_anal_tmax10_dt00001.npy")

# data1 = np.load("s3_r_gausss1_anal_tmax10_dt00001.npy")
# data2 = np.load("s3_r_gausss2_anal_tmax10_dt00001.npy")
# data3 = np.load("s3_r_gausss3_anal_tmax10_dt00001.npy")

# data1 = np.load("h3_r_gausss1_sim_tmax10_dt00001.npy")
# data2 = np.load("h3_r_gausss2_sim_tmax10_dt00001.npy")
# data3 = np.load("h3_r_gausss3_sim_tmax10_dt00001.npy")

data1 = np.load("s3_r_gausss1_sim_tmax10_dt1.npy")
data2 = np.load("s3_r_gausss2_sim_tmax10_dt1.npy")
data3 = np.load("s3_r_gausss3_sim_tmax10_dt1.npy")

# data1 = np.load("s3_r_gausss3_gt_tmax10_dt00001.npy")


fig,ax=plt.subplots(1,1)

distdata1 = np.zeros(np.shape(data1)[0])
distdata2 = np.zeros(np.shape(data2)[0])
distdata3 = np.zeros(np.shape(data3)[0])

counter = 0
for a in range(np.shape(data1)[0]):
    distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
    distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
    distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
    counter += 1

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
ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
ax.legend()
ax.set_title('Simulation Data')
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







