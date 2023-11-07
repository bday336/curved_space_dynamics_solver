# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from bumpgeofunc import bumpgeofunc
from bumpgeojac import bumpgeojac

def bumpgeodist(state_vec):

    def ppqfunc(a1, b1, g1, a2, b2, g2):
        return np.sqrt((a1 + a2)**2. + (b1 + b2)**2. + (g1 + g2)**2.)

    def pnqfunc(a1, b1, g1, a2, b2, g2):
        return np.sqrt((a1 - a2)**2. + (b1 - b2)**2. + (g1 - g2)**2.)

    a1,b1,g1,a2,b2,g2,ad1,bd1,gd1,ad2,bd2,gd2 = state_vec

    # PPQ & PNQ Functions
    ppq = ppqfunc(a1, b1, g1, a2, b2, g2)
    pnq = pnqfunc(a1, b1, g1, a2, b2, g2)

    fp = np.arcsinh(1/np.sqrt(2)*(ppq + pnq)/2) + (ppq + pnq)/2*np.sqrt(((ppq + pnq)/2)**2 + 2)/2
    fn = np.arcsinh(1/np.sqrt(2)*(ppq - pnq)/2) + (ppq - pnq)/2*np.sqrt(((ppq - pnq)/2)**2 + 2)/2

    return (fp - fn)

# Solver Setup

# Time array based on time step
dt = .01    # Number of steps
t_max = 50      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Time array based on number of steps
# nump = 10000    # Number of steps
# t_max = 10      # Total simulation time
# t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

# Simulation data container

gs1simdatalist = np.zeros((t_arr.shape[0],12))
gs2simdatalist = np.zeros((t_arr.shape[0],12))
gs3simdatalist = np.zeros((t_arr.shape[0],12))

startvec = np.array([
    [-25,.5,0.],[-25,-.5,0.],
    [1, (np.sqrt(2 + (25 + .5)**2) - np.sqrt(2 + (25 - .5)**2))/(np.sqrt(2 + (25 + .5)**2) + np.sqrt(2 + (25 - .5)**2)), 0], 
    [1, -(np.sqrt(2 + (25 + .5)**2) - np.sqrt(2 + (25 - .5)**2))/(np.sqrt(2 + (25 + .5)**2) + np.sqrt(2 + (25 - .5)**2)), 0]]).flatten()

# Initial Data
ks = 1.     # Spring Stiffness
x = bumpgeodist(startvec)      # Spring Rest Length
m = 1.      # Mass of point masses
params = [ks,x,m,m]

gs1simdatalist[0] = startvec.copy()
gs2simdatalist[0] = startvec.copy()
gs3simdatalist[0] = startvec.copy()

# First Step
step = 1

# Sim 
gs1simdatalist[step] = gausss1(startvec=startvec,params=params,dynfunc=bumpgeofunc,dynjac=bumpgeojac,dt=dt,tol=1e-15)
gs2simdatalist[step] = gausss2(startvec=startvec,params=params,dynfunc=bumpgeofunc,dynjac=bumpgeojac,dt=dt,tol=1e-15)
gs3simdatalist[step] = gausss3(startvec=startvec,params=params,dynfunc=bumpgeofunc,dynjac=bumpgeojac,dt=dt,tol=1e-15)

# Sim 
startvecgs1sim = gs1simdatalist[step]
startvecgs2sim = gs2simdatalist[step]
startvecgs3sim = gs3simdatalist[step]

step += 1

while (step <= int(t_max/dt)):

    # Sim 
    gs1simdatalist[step] = gausss1(startvec=startvecgs1sim,params=params,dynfunc=bumpgeofunc,dynjac=bumpgeojac,dt=dt,tol=1e-15)
    gs2simdatalist[step] = gausss2(startvec=startvecgs2sim,params=params,dynfunc=bumpgeofunc,dynjac=bumpgeojac,dt=dt,tol=1e-15)
    gs3simdatalist[step] = gausss3(startvec=startvecgs3sim,params=params,dynfunc=bumpgeofunc,dynjac=bumpgeojac,dt=dt,tol=1e-15)

    # Sim 
    startvecgs1sim = gs1simdatalist[step]
    startvecgs2sim = gs2simdatalist[step]
    startvecgs3sim = gs3simdatalist[step]

    if step%int(1/dt)==0:
            print(step)
    step += 1


# Sim in S3
np.save("bumpgeo_gausss1_sim_tmax10_dtp01",gs1simdatalist)
np.save("bumpgeo_gausss2_sim_tmax10_dtp01",gs2simdatalist)
np.save("bumpgeo_gausss3_sim_tmax10_dtp01",gs3simdatalist)


data1 = np.load("bumpgeo_gausss1_sim_tmax10_dtp01.npy")
data2 = np.load("bumpgeo_gausss2_sim_tmax10_dtp01.npy")
data3 = np.load("bumpgeo_gausss3_sim_tmax10_dtp01.npy")



fig,ax=plt.subplots(1,1)

distdata1 = np.zeros(np.shape(data1)[0])
distdata2 = np.zeros(np.shape(data2)[0])
distdata3 = np.zeros(np.shape(data3)[0])

counter = 0
for a in range(np.shape(data1)[0]):
    distdata1[counter] = bumpgeodist(data1[a])
    distdata2[counter] = bumpgeodist(data2[a])
    distdata3[counter] = bumpgeodist(data3[a])
    counter += 1


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







