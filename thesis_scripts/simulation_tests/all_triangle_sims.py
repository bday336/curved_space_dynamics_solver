# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3exacttriangle, dynjac_h3exacttriangle, dynfunc_h3simtriangle, dynjac_h3simtriangle, dynfunc_s3exacttriangle, dynjac_s3exacttriangle, dynfunc_s3simtriangle, dynjac_s3simtriangle

# Solver Setup

# Time array based on time step
dt = .00001    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Time array based on number of steps
# nump = 10000    # Number of steps
# t_max = 10      # Total simulation time
# t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

# Simulation data container

# Exact H3
# gs1exactdatalist = np.zeros((t_arr.shape[0],4))
# gs2exactdatalist = np.zeros((t_arr.shape[0],4))
# gs3exactdatalist = np.zeros((t_arr.shape[0],4))

# Sim H3
# gs1simdatalist = np.zeros((t_arr.shape[0],18))
# gs2simdatalist = np.zeros((t_arr.shape[0],18))
# gs3simdatalist = np.zeros((t_arr.shape[0],18))

# Exact S3
# gs1exactdatalist = np.zeros((t_arr.shape[0],4))
# gs2exactdatalist = np.zeros((t_arr.shape[0],4))
# gs3exactdatalist = np.zeros((t_arr.shape[0],4))

# Sim S3
gs1simdatalist = np.zeros((t_arr.shape[0],18))
gs2simdatalist = np.zeros((t_arr.shape[0],18))
gs3simdatalist = np.zeros((t_arr.shape[0],18))

# Initial Data
v = 1.      # Initial Velocity
ks = 1.     # Spring Stiffness
# x = 1.      # Spring Rest Length H3
# x = np.pi - 1.      # Spring Rest Length S3 exact
x = 1.      # Spring Rest Length S3
m = 1.      # Mass of point masses
params = [v,ks,x,m]

# Exact triangle in H3
# startvec = np.array([x/2.,np.arccosh(np.cosh(x)/np.cosh(x/2.)),0.,0.])
# Sim triangle in H3
# startvec = np.array([
#     [np.arccosh(np.cosh(x)/np.cosh(x/2.)),np.pi/2.,0.],[.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([np.arccosh(np.cosh(x)/np.cosh(x/2.)),np.pi/2.,0.],v,"x"), killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()
# Exact triangle in S3
# startvec = np.array([x/2.,np.arccos(np.cos(x)/np.cos(x/2.)),0.,0.])
# Sim bar in S3
startvec = np.array([
    [np.pi/2.,np.pi/2.,np.arccos(np.cos(x)/np.cos(x/2.))],[(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0],
    killingvecs3([np.pi/2.,np.pi/2.,np.arccos(np.cos(x)/np.cos(x/2.))],-v,"vz"), killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0],-v,"vz")]).flatten()

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

# Sim in H3
gs1simdatalist[0] = startvec.copy()
gs2simdatalist[0] = startvec.copy()
gs3simdatalist[0] = startvec.copy()

# First Step
step = 1

# Exact in H3
# gs1exactdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3exacttriangle,dynjac=dynjac_h3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
# gs2exactdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3exacttriangle,dynjac=dynjac_h3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
# gs3exactdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3exacttriangle,dynjac=dynjac_h3exacttriangle,dt=dt,tol=1e-15,imax=1e6)

# Sim in H3
# gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3simtriangle,dynjac=dynjac_h3simtriangle,dt=dt,tol=1e-15)
# gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3simtriangle,dynjac=dynjac_h3simtriangle,dt=dt,tol=1e-15)
# gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3simtriangle,dynjac=dynjac_h3simtriangle,dt=dt,tol=1e-15) 

# Exact in S3
# gs1exactdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3exacttriangle,dynjac=dynjac_s3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
# gs2exactdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3exacttriangle,dynjac=dynjac_s3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
# gs3exactdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3exacttriangle,dynjac=dynjac_s3exacttriangle,dt=dt,tol=1e-15,imax=1e6)

# Sim in S3
gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3simtriangle,dynjac=dynjac_s3simtriangle,dt=dt,tol=1e-15)
gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3simtriangle,dynjac=dynjac_s3simtriangle,dt=dt,tol=1e-15)
gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3simtriangle,dynjac=dynjac_s3simtriangle,dt=dt,tol=1e-15)

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

while (step <= int(t_max/dt+1)):
    # Exact in H3
    # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_h3exacttriangle,dynjac=dynjac_h3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
    # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_h3exacttriangle,dynjac=dynjac_h3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
    # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_h3exacttriangle,dynjac=dynjac_h3exacttriangle,dt=dt,tol=1e-15,imax=1e6)

    #Sim in H3
    # gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_h3simtriangle,dynjac=dynjac_h3simtriangle,dt=dt,tol=1e-15)
    # gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_h3simtriangle,dynjac=dynjac_h3simtriangle,dt=dt,tol=1e-15)
    # gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_h3simtriangle,dynjac=dynjac_h3simtriangle,dt=dt,tol=1e-15)    

    # Exact in S3
    # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_s3exacttriangle,dynjac=dynjac_s3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
    # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_s3exacttriangle,dynjac=dynjac_s3exacttriangle,dt=dt,tol=1e-15,imax=1e6)
    # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_s3exacttriangle,dynjac=dynjac_s3exacttriangle,dt=dt,tol=1e-15,imax=1e6)

    # Sim in S3
    gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_s3simtriangle,dynjac=dynjac_s3simtriangle,dt=dt,tol=1e-15)
    gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_s3simtriangle,dynjac=dynjac_s3simtriangle,dt=dt,tol=1e-15)
    gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_s3simtriangle,dynjac=dynjac_s3simtriangle,dt=dt,tol=1e-15)

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

    if step%int(1/dt+1)==0:
            print(step)
    step += 1

# Exact in H3
# np.save("h3_t_gausss1_anal_tmax10_dt00001",gs1exactdatalist)
# np.save("h3_t_gausss2_anal_tmax10_dt00001",gs2exactdatalist)
# np.save("h3_t_gausss3_anal_tmax10_dt00001",gs3exactdatalist)

# np.save("h3_t_gausss3_gt_tmax10_dt000005",gs3exactdatalist)

# Sim in H3
# np.save("h3_t_gausss1_sim_tmax10_dt1",gs1simdatalist)
# np.save("h3_t_gausss2_sim_tmax10_dt1",gs2simdatalist)
# np.save("h3_t_gausss3_sim_tmax10_dt1",gs3simdatalist)

# Exact in S3
# np.save("s3_t_gausss1_anal_tmax10_dt00001",gs1exactdatalist)
# np.save("s3_t_gausss2_anal_tmax10_dt00001",gs2exactdatalist)
# np.save("s3_t_gausss3_anal_tmax10_dt00001",gs3exactdatalist)

# np.save("s3_t_gausss3_gt_tmax10_dt00001",gs3exactdatalist)

# Sim in S3
np.save("s3_t_gausss1_sim_tmax10_dt00001",gs1simdatalist)
np.save("s3_t_gausss2_sim_tmax10_dt00001",gs2simdatalist)
np.save("s3_t_gausss3_sim_tmax10_dt00001",gs3simdatalist)

# data1 = np.load("h3_t_gausss1_anal_tmax10_dt00001.npy")
# data2 = np.load("h3_t_gausss2_anal_tmax10_dt00001.npy")
# data3 = np.load("h3_t_gausss3_anal_tmax10_dt00001.npy")

# data1 = np.load("s3_t_gausss1_anal_tmax10_dt00001.npy")
# data2 = np.load("s3_t_gausss2_anal_tmax10_dt00001.npy")
# data3 = np.load("s3_t_gausss3_anal_tmax10_dt00001.npy")

# data1 = np.load("h3_t_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("h3_t_gausss2_sim_tmax10_dt1.npy")
# data3 = np.load("h3_t_gausss3_sim_tmax10_dt1.npy")

data1 = np.load("s3_t_gausss1_sim_tmax10_dt00001.npy")
data2 = np.load("s3_t_gausss2_sim_tmax10_dt00001.npy")
data3 = np.load("s3_t_gausss3_sim_tmax10_dt00001.npy")

# data3 = np.load("s3_t_gausss3_gt_tmax10_dt00001.npy")

fig,ax=plt.subplots(1,1)

distdatasp12g1 = np.zeros(np.shape(data1)[0])
distdatasp12g2 = np.zeros(np.shape(data2)[0])
distdatasp12g3 = np.zeros(np.shape(data3)[0])

distdatasp13g1 = np.zeros(np.shape(data1)[0])
distdatasp13g2 = np.zeros(np.shape(data2)[0])
distdatasp13g3 = np.zeros(np.shape(data3)[0])

distdatasp23g1 = np.zeros(np.shape(data1)[0])
distdatasp23g2 = np.zeros(np.shape(data2)[0])
distdatasp23g3 = np.zeros(np.shape(data3)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
#     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
#     distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
#     counter += 1

counter = 0
for a in range(np.shape(data1)[0]):
    distdatasp12g1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
    distdatasp12g2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
    distdatasp12g3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))

    distdatasp13g1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][6:9]))
    distdatasp13g2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][6:9]))
    distdatasp13g3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][6:9]))

    distdatasp23g1[counter] = r4dist(rot2r4(data1[a][3:6]),rot2r4(data1[a][6:9]))
    distdatasp23g2[counter] = r4dist(rot2r4(data2[a][3:6]),rot2r4(data2[a][6:9]))
    distdatasp23g3[counter] = r4dist(rot2r4(data3[a][3:6]),rot2r4(data3[a][6:9]))
    counter += 1

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdatasp12g1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
#     distdatasp12g2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
#     distdatasp12g3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))

#     distdatasp13g1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][6:9]))
#     distdatasp13g2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][6:9]))
#     distdatasp13g3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][6:9]))

#     distdatasp23g1[counter] = h3dist(rot2hyp(data1[a][3:6]),rot2hyp(data1[a][6:9]))
#     distdatasp23g2[counter] = h3dist(rot2hyp(data2[a][3:6]),rot2hyp(data2[a][6:9]))
#     distdatasp23g3[counter] = h3dist(rot2hyp(data3[a][3:6]),rot2hyp(data3[a][6:9]))
#     counter += 1

# ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")

ax.plot(t_arr,distdatasp12g1,'r',label = "Gauss s1")
ax.plot(t_arr,distdatasp12g2,'r',label = "Gauss s2")
ax.plot(t_arr,distdatasp12g3,'r',label = "Gauss s3")

ax.plot(t_arr,distdatasp13g1,'k',label = "Gauss s1")
ax.plot(t_arr,distdatasp13g2,'k',label = "Gauss s2")
ax.plot(t_arr,distdatasp13g3,'k',label = "Gauss s3")

ax.plot(t_arr,distdatasp23g1,'b',label = "Gauss s1")
ax.plot(t_arr,distdatasp23g2,'b',label = "Gauss s2")
ax.plot(t_arr,distdatasp23g3,'b',label = "Gauss s3")

# ax.plot(t_arr,2*data3[:,0],'b',label = "Gauss s3")
# ax.plot(t_arr,np.arccosh(np.cosh(data3[:,0])*np.cosh(data3[:,1])),'r',label = "Gauss s3")

# ax.plot(t_arr,2*data3[:,0],'b',label = "Gauss s3")
# ax.plot(t_arr,np.arccos(np.cos(data3[:,0])*np.cos(data3[:,1])),'r',label = "Gauss s3")
ax.legend()
ax.set_title('Simulation Data')
ax.set_ylim(0,2.5)
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







