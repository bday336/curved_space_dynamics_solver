# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from derivative_bank import *
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3, rads2dae
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
# gs1simdatalist = np.zeros((t_arr.shape[0],12+1))
rs2simdatalist = np.zeros((t_arr.shape[0],12+1))
# gs3simdatalist = np.zeros((t_arr.shape[0],12+1))

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
# x = 1.      # Spring Rest Length H3
x = 1.      # Rod Length H3
m1 = 1.      # Mass of point masses
m2 = 1.      # Mass of point masses
params = [m1,m2,x]

# Exact bar in H3
# startvec = np.array([.5,0.])
# Sim bar in H3
startvec = np.array([
    [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
    killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()
startvec = np.append(startvec,0.)
# Exact bar in S3
# startvec = np.array([(np.pi - 1.)/2.,0.])
# Sim bar in S3
# startvec = np.array([
#     [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
#     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()



# Exact in H3
# gs1exactdatalist[0] = startvec.copy()
# gs2exactdatalist[0] = startvec.copy()
# gs3exactdatalist[0] = startvec.copy()

# Sim in H3
# gs1simdatalist[0] = startvec.copy()
rs2simdatalist[0] = startvec.copy()
# gs3simdatalist[0] = startvec.copy()

# Exact in S3
# gs1exactdatalist[0] = startvec.copy()
# gs2exactdatalist[0] = startvec.copy()
# gs3exactdatalist[0] = startvec.copy()

# Sim in H3
# gs1simdatalist[0] = startvec.copy()
# gs2simdatalist[0] = startvec.copy()
# gs3simdatalist[0] = startvec.copy()

# First Step
step = 1

# stepn = startvec.copy()
# stepc1 = startvec.copy()
# stepc2 = startvec.copy()
# initvec = np.array([
#     stepc1[0],stepc2[0],
#     stepc1[1],stepc2[1],
#     stepc1[2],stepc2[2],

#     stepc1[3],stepc2[3],
#     stepc1[4],stepc2[4],
#     stepc1[5],stepc2[5],

#     stepc1[6],stepc2[6],
#     stepc1[7],stepc2[7],
#     stepc1[8],stepc2[8],

#     stepc1[9],stepc2[9],
#     stepc1[10],stepc2[10],
#     stepc1[11],stepc2[11],

#     stepc1[12],stepc2[12]
#     ])

# ex1,dex1 = h3rads2rodex1(stepn,stepc1,stepc2,dt,params)
# ex2,dex2 = h3rads2rodex2(stepn,stepc1,stepc2,dt,params)
# ex3,dex3 = h3rads2rodex3(stepn,stepc1,stepc2,dt,params)
# ex4,dex4 = h3rads2rodex4(stepn,stepc1,stepc2,dt,params)
# ex5,dex5 = h3rads2rodex5(stepn,stepc1,stepc2,dt,params)
# ex6,dex6 = h3rads2rodex6(stepn,stepc1,stepc2,dt,params)

# ex7,dex7 = h3rads2rodex7(stepn,stepc1,stepc2,dt,params)
# ex8,dex8 = h3rads2rodex8(stepn,stepc1,stepc2,dt,params)
# ex9,dex9 = h3rads2rodex9(stepn,stepc1,stepc2,dt,params)
# ex10,dex10 = h3rads2rodex10(stepn,stepc1,stepc2,dt,params)
# ex11,dex11 = h3rads2rodex11(stepn,stepc1,stepc2,dt,params)
# ex12,dex12 = h3rads2rodex12(stepn,stepc1,stepc2,dt,params)

# ex13,dex13 = h3rads2rodex13(stepn,stepc1,stepc2,dt,params)
# ex14,dex14 = h3rads2rodex14(stepn,stepc1,stepc2,dt,params)
# ex15,dex15 = h3rads2rodex15(stepn,stepc1,stepc2,dt,params)
# ex16,dex16 = h3rads2rodex16(stepn,stepc1,stepc2,dt,params)
# ex17,dex17 = h3rads2rodex17(stepn,stepc1,stepc2,dt,params)
# ex18,dex18 = h3rads2rodex18(stepn,stepc1,stepc2,dt,params)

# ex19,dex19 = h3rads2rodex19(stepn,stepc1,stepc2,dt,params)
# ex20,dex20 = h3rads2rodex20(stepn,stepc1,stepc2,dt,params)
# ex21,dex21 = h3rads2rodex21(stepn,stepc1,stepc2,dt,params)
# ex22,dex22 = h3rads2rodex22(stepn,stepc1,stepc2,dt,params)
# ex23,dex23 = h3rads2rodex23(stepn,stepc1,stepc2,dt,params)
# ex24,dex24 = h3rads2rodex24(stepn,stepc1,stepc2,dt,params)

# ex25,dex25 = h3rads2rodex25(stepn,stepc1,stepc2,dt,params)
# ex26,dex26 = h3rads2rodex26(stepn,stepc1,stepc2,dt,params)

# jacobian = np.array([
#     dex1,dex2,dex3,
#     dex4,dex5,dex6,

#     dex7,dex8,dex9,
#     dex10,dex11,dex12,

#     dex13,dex14,dex15,
#     dex16,dex17,dex18,

#     dex19,dex20,dex21,
#     dex22,dex23,dex24,

#     dex25,dex26

# ])

# conlist = np.array([
#     ex1,ex2,ex3,
#     ex4,ex5,ex6,

#     ex7,ex8,ex9,
#     ex10,ex11,ex12,

#     ex13,ex14,ex15,
#     ex16,ex17,ex18,

#     ex19,ex20,ex21,
#     ex22,ex23,ex24,

#     ex25,ex26
# ])

# diff1 = np.linalg.solve(jacobian,-conlist)

# val1 = diff1 + initvec

# # Begin Iterations
# counter = 0
# while (np.linalg.norm(diff1) >= 1e-10 and counter <= 100):
#     stepc1 = val1[0::2].copy()
#     stepc2 = val1[1::2].copy()

#     ex1,dex1 = h3rads2rodex1(stepn,stepc1,stepc2,dt,params)
#     ex2,dex2 = h3rads2rodex2(stepn,stepc1,stepc2,dt,params)
#     ex3,dex3 = h3rads2rodex3(stepn,stepc1,stepc2,dt,params)
#     ex4,dex4 = h3rads2rodex4(stepn,stepc1,stepc2,dt,params)
#     ex5,dex5 = h3rads2rodex5(stepn,stepc1,stepc2,dt,params)
#     ex6,dex6 = h3rads2rodex6(stepn,stepc1,stepc2,dt,params)

#     ex7,dex7 = h3rads2rodex7(stepn,stepc1,stepc2,dt,params)
#     ex8,dex8 = h3rads2rodex8(stepn,stepc1,stepc2,dt,params)
#     ex9,dex9 = h3rads2rodex9(stepn,stepc1,stepc2,dt,params)
#     ex10,dex10 = h3rads2rodex10(stepn,stepc1,stepc2,dt,params)
#     ex11,dex11 = h3rads2rodex11(stepn,stepc1,stepc2,dt,params)
#     ex12,dex12 = h3rads2rodex12(stepn,stepc1,stepc2,dt,params)

#     ex13,dex13 = h3rads2rodex13(stepn,stepc1,stepc2,dt,params)
#     ex14,dex14 = h3rads2rodex14(stepn,stepc1,stepc2,dt,params)
#     ex15,dex15 = h3rads2rodex15(stepn,stepc1,stepc2,dt,params)
#     ex16,dex16 = h3rads2rodex16(stepn,stepc1,stepc2,dt,params)
#     ex17,dex17 = h3rads2rodex17(stepn,stepc1,stepc2,dt,params)
#     ex18,dex18 = h3rads2rodex18(stepn,stepc1,stepc2,dt,params)

#     ex19,dex19 = h3rads2rodex19(stepn,stepc1,stepc2,dt,params)
#     ex20,dex20 = h3rads2rodex20(stepn,stepc1,stepc2,dt,params)
#     ex21,dex21 = h3rads2rodex21(stepn,stepc1,stepc2,dt,params)
#     ex22,dex22 = h3rads2rodex22(stepn,stepc1,stepc2,dt,params)
#     ex23,dex23 = h3rads2rodex23(stepn,stepc1,stepc2,dt,params)
#     ex24,dex24 = h3rads2rodex24(stepn,stepc1,stepc2,dt,params)

#     ex25,dex25 = h3rads2rodex25(stepn,stepc1,stepc2,dt,params)
#     ex26,dex26 = h3rads2rodex26(stepn,stepc1,stepc2,dt,params)

#     jacobian = np.array([
#         dex1,dex2,dex3,
#         dex4,dex5,dex6,

#         dex7,dex8,dex9,
#         dex10,dex11,dex12,

#         dex13,dex14,dex15,
#         dex16,dex17,dex18,

#         dex19,dex20,dex21,
#         dex22,dex23,dex24,

#         dex25,dex26

#     ])

#     conlist = np.array([
#         ex1,ex2,ex3,
#         ex4,ex5,ex6,

#         ex7,ex8,ex9,
#         ex10,ex11,ex12,

#         ex13,ex14,ex15,
#         ex16,ex17,ex18,

#         ex19,ex20,ex21,
#         ex22,ex23,ex24,

#         ex25,ex26
#     ])

#     diff2 = np.linalg.solve(jacobian,-conlist)

#     val2 = diff2 + val1

#     val1 = val2
#     diff1 = diff2
#     counter += 1

# Exact in H3
# gs1exactdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)
# gs2exactdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)
# gs3exactdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)

# Sim in H3
# gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt)
rs2simdatalist[step] = rads2dae(startvec=startvec,params=params,dt=dt)
# gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt) 

# Exact in S3
# gs1exactdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)
# gs2exactdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)
# gs3exactdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)

# Sim in S3
# gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt)
# gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt)
# gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt) 

# Exact in H3
# startvecgs1ex = gs1exactdatalist[step]
# startvecgs2ex = gs2exactdatalist[step]
# startvecgs3ex = gs3exactdatalist[step]

# Sim in H3
# startvecgs1sim = gs1simdatalist[step]
startvecrs2sim = rs2simdatalist[step]
# startvecgs3sim = gs3simdatalist[step]

# Exact in S3
# startvecgs1ex = gs1exactdatalist[step]
# startvecgs2ex = gs2exactdatalist[step]
# startvecgs3ex = gs3exactdatalist[step]

# Sim in S3
# startvecgs1sim = gs1simdatalist[step]
# startvecgs2sim = gs2simdatalist[step]
# startvecgs3sim = gs3simdatalist[step]

step += 1

while (step <= int(t_max/dt)):
    # Exact in H3
    # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)
    # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)
    # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)

    # Sim in H3
    # gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt)
    rs2simdatalist[step] = rads2dae(startvec=startvecrs2sim,params=params,dt=dt)
    # gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt)    

    # Exact in S3
    # gs1exactdatalist[step] = gausss1(startvec=startvecgs1ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)
    # gs2exactdatalist[step] = gausss2(startvec=startvecgs2ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)
    # gs3exactdatalist[step] = gausss3(startvec=startvecgs3ex,params=params,dynfunc=dynfunc_s3exactbar,dynjac=dynjac_s3exactbar,dt=dt)

    # Sim in S3
    # gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt)
    # gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt)
    # gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=dynfunc_s3simbar,dynjac=dynjac_s3simbar,dt=dt) 

    # Exact in H3
    # startvecgs1ex = gs1exactdatalist[step]
    # startvecgs2ex = gs2exactdatalist[step]
    # startvecgs3ex = gs3exactdatalist[step]

    # Sim in H3
    # startvecgs1sim = gs1simdatalist[step]
    startvecrs2sim = rs2simdatalist[step]
    # startvecgs3sim = gs3simdatalist[step]

    # Exact in S3
    # startvecgs1ex = gs1exactdatalist[step]
    # startvecgs2ex = gs2exactdatalist[step]
    # startvecgs3ex = gs3exactdatalist[step]

    # Sim in S3
    # startvecgs1sim = gs1simdatalist[step]
    # startvecgs2sim = gs2simdatalist[step]
    # startvecgs3sim = gs3simdatalist[step]

    if step%int(1/dt)==0:
            print(step)
    step += 1

# Exact in H3
# np.save("gausss1_tmax10_dt01",gs1exactdatalist)
# np.save("gausss2_tmax10_dt01",gs2exactdatalist)
# np.save("h3_r_gausss3_gt_tmax10_dt000005",gs3exactdatalist)

# Sim in H3
# np.save("h3_r_gausss1_sim_tmax10_dt00001",gs1simdatalist)
# np.save("h3_r_gausss2_sim_tmax10_dt00001",gs2simdatalist)
# np.save("h3_r_gausss3_sim_tmax10_dt00001",gs3simdatalist)

# Exact in S3
# np.save("gausss1_tmax10_dt01",gs1exactdatalist)
# np.save("gausss2_tmax10_dt01",gs2exactdatalist)
# np.save("gausss3_gt_tmax10_dt000005",gs3exactdatalist)

# Sim in S3
# np.save("s3_gausss1_sim_tmax10_dt00001",gs1simdatalist)
# np.save("s3_gausss2_sim_tmax10_dt00001",gs2simdatalist)
# np.save("s3_gausss3_sim_tmax10_dt00001",gs3simdatalist)

# data1 = np.load("h3_r_gausss1_sim_tmax10_dt00001.npy")
# data2 = np.load("h3_r_gausss2_sim_tmax10_dt00001.npy")
# data3 = np.load("h3_r_gausss3_sim_tmax10_dt00001.npy")

# data1 = np.load("h3_r_gausss3_gt_tmax10_dt000005.npy")

data2 = rs2simdatalist

fig,ax=plt.subplots(1,1)

# distdata1 = np.zeros(np.shape(data1)[0])
distdata2 = np.zeros(np.shape(data2)[0])
# distdata3 = np.zeros(np.shape(data3)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
#     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
#     distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
#     counter += 1

counter = 0
for a in range(np.shape(data2)[0]):
    # distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
    distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
    # distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
    counter += 1

# ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")
# ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
ax.legend()
ax.set_title('Simulation Data')
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







