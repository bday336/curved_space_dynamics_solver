# Allow for package importing
import os, sys
sys.path.append(os.getcwd()) 

import numpy as np
import matplotlib.pyplot as plt
from build.src.function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from build.src.integrator_files.integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from build.src.test_system_simulations.test_system_bank import dynfunc_h3exactbar, dynjac_h3exactbar, dynfunc_h3simbar, dynjac_h3simbar, dynfunc_s3exactbar, dynjac_s3exactbar, dynfunc_s3simbar, dynjac_s3simbar



# Initialize Test System
                
dt = .1         # Number of steps
t_max = 10.      # Total simulation time

# Initial Data
v = 1.      # Initial Velocity
ks = 1.     # Spring Stiffness
x = 1.      # Spring Rest Length H3 and S3
# x = np.pi - 1.      # Spring Rest Length S3 exact
m = 1.      # Mass of point masses
params = [v,ks,x,m]


# Exact bar in H3
# geometry, system_id = "h3", "exact"
# startvec = np.array([.5,0.])
# Sim bar in H3
geometry = "h3"
startvec = np.array([
    [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
    killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()

# Exact bar in S3
# geometry, system_id = "s3", "exact"
# startvec = np.array([(np.pi - 1.)/2.,0.])
# Sim bar in S3
# geometry, system_id = "s3", "sim"
# startvec = np.array([
#     [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
#     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()

# Solver 
solver_id = "gs1"

# Initialize Simulation Object
sim_test = SpringBarSimulation(geometry, params, dt, t_max, solver_id, True)
sim_test.set_initial_conditions(startvec)
sim_test.run()

data1 = np.load("h3_r_gs1_sim_tmax10.0_dt0.1.npy")


# data1 = np.load("h3_r_gausss1_anal_tmax10_dt00001.npy")
# data2 = np.load("h3_r_gausss2_anal_tmax10_dt00001.npy")
# data3 = np.load("h3_r_gausss3_anal_tmax10_dt00001.npy")

# data1 = np.load("s3_r_gausss1_anal_tmax10_dt00001.npy")
# data2 = np.load("s3_r_gausss2_anal_tmax10_dt00001.npy")
# data3 = np.load("s3_r_gausss3_anal_tmax10_dt00001.npy")

# data1 = np.load("h3_r_gausss1_sim_tmax10_dt00001.npy")
# data2 = np.load("h3_r_gausss2_sim_tmax10_dt00001.npy")
# data3 = np.load("h3_r_gausss3_sim_tmax10_dt00001.npy")

# data1 = np.load("s3_r_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("s3_r_gausss2_sim_tmax10_dt1.npy")
# data3 = np.load("s3_r_gausss3_sim_tmax10_dt1.npy")

# # data1 = np.load("s3_r_gausss3_gt_tmax10_dt00001.npy")


fig,ax=plt.subplots(1,1)

distdata1 = np.zeros(np.shape(data1)[0])
# distdata2 = np.zeros(np.shape(data2)[0])
# distdata3 = np.zeros(np.shape(data3)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
#     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
#     distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
#     counter += 1

counter = 0
for a in range(np.shape(data1)[0]):
    distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
    # distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
    # distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
    counter += 1

# ax.plot(t_arr,2.*(np.pi/2. - data1[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - data2[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - data3[:,0]),'b',label = "Gauss s3")
# ax.plot(sim_test.t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,2.*(data2[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,2.*(data3[:,0]),'b',label = "Gauss h3")
ax.plot(sim_test.t_arr,distdata1,'r',label = "Gauss s1")
# ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
ax.legend()
ax.set_title('Simulation Data')
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()








        
    