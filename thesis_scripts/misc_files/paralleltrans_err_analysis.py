# Generic Simulation Script

import numpy as np
from numpy import sinh,cosh,sin,cos,tanh,tan
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from test_system_bank import dynfunc_h3simgeo, dynjac_h3simgeo, dynfunc_s3simgeo, dynjac_s3simgeo

# Two point geodesic functions
def h32pgeotangent(p1,p2,dt):
    x1,y1,z1,w1 = rot2hyp(p1)
    x2,y2,z2,w2 = rot2hyp(p2)

    Q = x1*x2 + y1*y2 + z1*z2 - w1*w2

    # normat1 = np.sqrt(abs( x2+Q*x1 + y2+Q*y1 + z2+Q*z1 - w2+Q*w1))

    geodir = np.array([x2+Q*x1 , y2+Q*y1 , z2+Q*z1 , w2+Q*w1])/dt #/normat1
    return geodir

def s32pgeotangent(p1,p2,dt):
    x1,y1,z1,w1 = rot2r4(p1)
    x2,y2,z2,w2 = rot2r4(p2)

    Q = x1*x2 + y1*y2 + z1*z2 + w1*w2

    # normat1 = np.sqrt( x2+Q*x1 + y2+Q*y1 + z2+Q*z1 + w2+Q*w1)

    geodir = np.array([x2-Q*x1 , y2-Q*y1 , z2-Q*z1 , w2-Q*w1])/dt #/normat1
    return geodir

# Velocity transformation functions between parameters and embedding
def rot2hypv(pos, vel): 
    return np.array([
        vel[0]*cosh(pos[0])*sin(pos[1])*cos(pos[2]) + vel[1]*sinh(pos[0])*cos(pos[1])*cos(pos[2]) - vel[2]*sinh(pos[0])*sin(pos[1])*sin(pos[2]),
        vel[0]*cosh(pos[0])*sin(pos[1])*sin(pos[2]) + vel[1]*sinh(pos[0])*cos(pos[1])*sin(pos[2]) + vel[2]*sinh(pos[0])*sin(pos[1])*cos(pos[2]),
        vel[0]*cosh(pos[0])*cos(pos[1]) - vel[1]*sinh(pos[0])*sin(pos[1]),
        vel[0]*sinh(pos[0])])

def hyp2rotv(pos, vel): 
    return np.array([
        vel[3]/np.sinh(pos[0]),
        (vel[3]/np.tanh(pos[0])*np.cos(pos[1]) - vel[2])/(np.sinh(pos[0])*np.sin(pos[1])),
        (vel[1]*np.cos(pos[2]) - vel[0]*np.sin(pos[2]))/(np.sinh(pos[0])*np.sin(pos[1]))])

def rot2r4v(pos, vel): 
    return np.array([
        vel[0]*cos(pos[0])*sin(pos[1])*cos(pos[2]) + vel[1]*sin(pos[0])*cos(pos[1])*cos(pos[2]) - vel[2]*sin(pos[0])*sin(pos[1])*sin(pos[2]),
        vel[0]*cos(pos[0])*sin(pos[1])*sin(pos[2]) + vel[1]*sin(pos[0])*cos(pos[1])*sin(pos[2]) + vel[2]*sin(pos[0])*sin(pos[1])*cos(pos[2]),
        vel[0]*cos(pos[0])*cos(pos[1]) - vel[1]*sin(pos[0])*sin(pos[1]),
        -vel[0]*sin(pos[0])])

def r42rotv(pos, vel): 
    return np.array([
        -vel[3]/np.sin(pos[0]),
        (vel[3]/np.tan(pos[0])*np.cos(pos[1]) - vel[2])/(np.sin(pos[0])*np.sin(pos[1])),
        (vel[1]*np.cos(pos[2]) - vel[0]*np.sin(pos[2]))/(np.sin(pos[0])*np.sin(pos[1]))])

# Discrete Parallel Transport Step
def h3disptransport(data1,data2,dt):
    a1,b1,g1,ad1,bd1,gd1 = data1
    a2,b2,g2,ad2,bd2,gd2 = data2
    geodir = h32pgeotangent([a1,b1,g1],[a2,b2,g2],dt)
    geotan = hyp2rotv([a1,b1,g1],geodir)
    v_init = np.array([ad1,bd1,gd1])
    christterma = sinh(a1)*cosh(a1)*geotan[1]*bd1 + sinh(a1)*cosh(a1)*sin(b1)**2.*geotan[2]*gd1
    christtermb = -1./tanh(a1)*geotan[1]*ad1 - 1./tanh(a1)*geotan[0]*bd1 + sin(b1)*cos(b1)*geotan[2]*gd1
    christtermg = -1./tanh(a1)*geotan[2]*ad1 - 1./tanh(a1)*geotan[0]*gd1 - 1./tan(b1)*geotan[2]*bd1 - 1./tan(b1)*geotan[1]*gd1
    v_transp = v_init + np.array([dt*christterma,dt*christtermb,dt*christtermg])
    simvhyp = rot2hypv(data1[:3],data2[3:6])
    appvhyp = rot2hypv(data1[:3],v_transp)
    return appvhyp - simvhyp

def s3disptransport(data1,data2,dt):
    a1,b1,g1,ad1,bd1,gd1 = data1
    a2,b2,g2,ad2,bd2,gd2 = data2
    geodir = s32pgeotangent([a1,b1,g1],[a2,b2,g2],dt)
    geotan = r42rotv([a1,b1,g1],geodir)
    v_init = np.array([ad1,bd1,gd1])
    christterma = sin(a1)*cos(a1)*geotan[1]*bd1 + sin(a1)*cos(a1)*sin(b1)**2.*geotan[2]*gd1
    christtermb = -1./tan(a1)*geotan[1]*ad1 - 1./tan(a1)*geotan[0]*bd1 + sin(b1)*cos(b1)*geotan[2]*gd1
    christtermg = -1./tan(a1)*geotan[2]*ad1 - 1./tan(a1)*geotan[0]*gd1 - 1./tan(b1)*geotan[2]*bd1 - 1./tan(b1)*geotan[1]*gd1
    v_transp = v_init + np.array([dt*christterma,dt*christtermb,dt*christtermg])
    simvr4 = rot2r4v(data1[:3],data2[3:6])
    appvr4 = rot2r4v(data1[:3],v_transp)
    return appvr4 - simvr4

# Solver Setup

# # Time array based on time step
dt = .00001    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max,dt)

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
# gs1simdatalist = np.zeros((t_arr.shape[0],6))
# gs2simdatalist = np.zeros((t_arr.shape[0],6))
# gs3simdatalist = np.zeros((t_arr.shape[0],6))

# Initial Data
v = 1.      # Initial Velocity
m = 1.      # Mass of point masses
params = [v,m]

# Sim bar in H3
# startvec = np.array([
#     [.5,np.pi/2.,np.pi/2.],
#     killingvech3([.5,np.pi/2.,np.pi/2.],v,"x")]).flatten()
# startvec = np.array([
#     [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()

# Sim bar in S3
startvec = np.array([
    [(np.pi - 1.)/2.,np.pi/2.,0.],
    killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()
# startvec = np.array([
#     [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
#     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()



# data1 = np.load("h3_geo_gausss1_sim_tmax10_dt01.npy")
# data2 = np.load("h3_geo_gausss2_sim_tmax10_dt01.npy")
# data3 = np.load("h3_geo_gausss3_sim_tmax10_dt01.npy")

data1 = np.load("s3_geo_gausss1_sim_tmax10_dt00001.npy")
data2 = np.load("s3_geo_gausss2_sim_tmax10_dt00001.npy")
data3 = np.load("s3_geo_gausss3_sim_tmax10_dt00001.npy")

# data1 = np.load("s3_r_gausss3_gt_tmax10_dt00001.npy")


fig,ax=plt.subplots(1,1)

distdata1 = np.zeros(np.shape(data1)[0])
distdata2 = np.zeros(np.shape(data2)[0])
distdata3 = np.zeros(np.shape(data3)[0])

counter = 0
for a in range(np.shape(data1)[0]-1):
    distdata1[counter] = np.linalg.norm(s3disptransport(data1[a],data1[a+1],dt))
    distdata2[counter] = np.linalg.norm(s3disptransport(data2[a],data2[a+1],dt))
    distdata3[counter] = np.linalg.norm(s3disptransport(data3[a],data3[a+1],dt))
    counter += 1

# counter = 0
# for a in range(np.shape(data1)[0]-1):
#     distdata1[counter] = np.linalg.norm(h3disptransport(data1[a],data1[a+1],dt))
#     distdata2[counter] = np.linalg.norm(h3disptransport(data2[a],data2[a+1],dt))
#     distdata3[counter] = np.linalg.norm(h3disptransport(data3[a],data3[a+1],dt))
#     counter += 1

# ax.plot(t_arr,2.*(np.pi/2. - data1[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - data2[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - data3[:,0]),'b',label = "Gauss s3")
# ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,2.*(data2[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,2.*(data3[:,0]),'b',label = "Gauss h3")
ax.plot(t_arr,distdata1[:counter-1],'r',label = "Gauss s1")
ax.plot(t_arr,distdata2[:counter-1],'k',label = "Gauss s2")
ax.plot(t_arr,distdata3[:counter-1],'b',label = "Gauss s3")
ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
plt.show()







