# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3exactbar, dynjac_h3exactbar, dynfunc_h3simbar, dynjac_h3simbar

# Solver Setup

# Time array based on time step
dt = .000005    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

data1 = np.load("gausss1_sim_tmax10_dt00001.npy")
data2 = np.load("gausss2_sim_tmax10_dt00001.npy")
data3 = np.load("gausss3_sim_tmax10_dt00001.npy")


# fig,ax=plt.subplots(1,1)

# distdata1 = np.zeros(np.shape(data1)[0])
# distdata2 = np.zeros(np.shape(data2)[0])
# distdata3 = np.zeros(np.shape(data3)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
#     distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
#     distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
#     counter += 1

# ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
# ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()