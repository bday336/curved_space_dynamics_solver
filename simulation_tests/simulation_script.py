# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3exactbar, dynjac_h3exactbar, dynfunc_h3simbar, dynjac_h3simbar

# Solver Setup
nump = 10000    # Number of steps
t_max = 10      # Total simulation time
t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

datalist = np.zeros((t_arr.shape[0],12))     # Simulation data container

# Initial Data
v = 1.      # Initial Velocity
ks = 1.     # Spring Stiffness
x = 1.      # Spring Rest Length
m = 1.      # Mass of point masses
params = [v,ks,x,m]

# startvec = np.array([.5,0.])
startvec = np.array([
    [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
    killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()

datalist[0] = startvec.copy()

# First Step
step = 1

while (step <= t_max/dt):
    datalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3simbar,dynjac=dynjac_h3simbar,dt=dt)
    # datalist[step] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3exactbar,dynjac=dynjac_h3exactbar,dt=dt)
    startvec = datalist[step]
    if step%100==0:
            print(step)
    step += 1


fig,ax=plt.subplots(1,1)

distdata = np.zeros(np.shape(datalist)[0])
counter = 0
for a in datalist:
    distdata[counter] = h3dist(rot2hyp(a[0:3]),rot2hyp(a[3:6]))
    counter += 1

ax.plot(t_arr,distdata,'r')
# ax.plot(t_arr,datalist[0:,0],'r')
ax.set_title('Simulation Data')
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







