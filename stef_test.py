# # -*- coding: utf-8 -*-
# """
# Spyder Editor

# This is a temporary script file.
# """
# import math
# import numpy as np
# import matplotlib.pyplot as plt
# from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
# from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
# from collision_function_bank import h3collision
# from test_system_bank import dynfunc_h3sim2ballcol, dynjac_h3sim2ballcol, dynfunc_s3simbar, dynjac_s3simbar

# def f(x):
#     return (x*x)-(2*x)-56

# def my_f(dt_boo):
#     precollision = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_boo)
#     distck3 = h3dist(rot2hyp(precollision[0:3]),rot2hyp(precollision[3:6]))
#     return (distck3 - (r1 +r2))


# p0 = 0.1/2.
# tol = 0.0005
# n0 = 10 #000

# # Initial Data
# v = .5      # Initial Velocity
# x = 1.      # Spring Rest Length H3/S3
# m1 = 1.      # Mass of point mass 1
# m2 = 1.      # Mass of point mass 2
# r1 = .5      # Radius of mass 1
# r2 = .5      # Radisu of mass 2
# params = [m1,m2,r1,r2]


# # Sim bar in H3
# # startvec = np.array([
# #     [.5,np.pi/2.,0.*np.pi/2.],[.5,np.pi/2.,2.*np.pi/2.],
# #     killingvech3([.5,np.pi/2.,0.*np.pi/2.],-v,"x"), killingvech3([.5,np.pi/2.,2.*np.pi/2.],v,"x")]).flatten()

# # Mathematica check
# startvec = np.array([
#     [2.,np.pi/2.,4.*np.pi/4.],[2.,np.pi/2.,0.*np.pi/4.],
#     killingvech3([2.,np.pi/2.,4.*np.pi/4.],v,"x"), killingvech3([2.,np.pi/2.,0.*np.pi/4.],-v,"x")]).flatten()

# def Steffensen(f,p0,tol):
   
#     for i in range(1,n0) :
#         p1 = p0 + f(p0)
#         p2 = p1 + f(p1)
#         print(pow((p2 - p1),2))
#         print(p2 - (2*p1) + p0)
#         p = p2 - (pow((p2 - p1),2)/(p2 - (2*p1) + p0))
#         print(p-p0)
#         if abs(p-p0) < tol:
#             print("Converge after %f iterations"%i)
#             return p
#         print(p)
#         p0 = p
#     print('failed to converge in %f iterations' %n0)
#     return p

# ans = Steffensen(my_f,p0,tol)
# print("value is: %f"%ans)

# print("Checking Result: %f"%f(-6.549831825767762))










# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from collision_function_bank import h3collision
from test_system_bank import dynfunc_h3sim2ballcol, dynjac_h3sim2ballcol, dynfunc_s3simbar, dynjac_s3simbar

# Solver Setup

# Time array based on time step
dt = .001    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)


# Simulation data container

gs3simdatalist = np.zeros((t_arr.shape[0],12))


# Initial Data
v = .5      # Initial Velocity
x = 1.      # Spring Rest Length H3/S3
m1 = 1.      # Mass of point mass 1
m2 = 1.      # Mass of point mass 2
r1 = .5      # Radius of mass 1
r2 = .5      # Radisu of mass 2
params = [m1,m2,r1,r2]


# Sim bar in H3
# startvec = np.array([
#     [.5,np.pi/2.,0.*np.pi/2.],[.5,np.pi/2.,2.*np.pi/2.],
#     killingvech3([.5,np.pi/2.,0.*np.pi/2.],-v,"x"), killingvech3([.5,np.pi/2.,2.*np.pi/2.],v,"x")]).flatten()

# Sim bar in H3
# startvec = np.array([
#     [1.5,np.pi/2.,0.*np.pi/2.],[1.5,np.pi/2.,2.*np.pi/2.],
#     killingvech3([1.5,np.pi/2.,0.*np.pi/2.],-v,"x"), killingvech3([1.5,np.pi/2.,2.*np.pi/2.],v,"x")]).flatten()

# Mathematica check
startvec = np.array([
    [2.,np.pi/2.,2.*np.pi/4.],[2.,np.pi/2.,0.*np.pi/4.],
    killingvech3([2.,np.pi/2.,2.*np.pi/4.],-v,"y"), killingvech3([2.,np.pi/2.,0.*np.pi/4.],-v,"x")]).flatten()


# Initial condition t=0 (step=0)

gs3simdatalist[0] = startvec.copy()
distck3 = h3dist(rot2hyp(gs3simdatalist[0][0:3]),rot2hyp(gs3simdatalist[0][3:6]))
print("step 0")
print(distck3)


# First Step
step = 1
gs3simdatalist[step] = gausss3(startvec=gs3simdatalist[step-1],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt) 
distck3 = h3dist(rot2hyp(gs3simdatalist[step][0:3]),rot2hyp(gs3simdatalist[step][3:6]))
print("step 1")
print(distck3)

# if distck3 < r1 + r2:
#     print("Collided")
#     tol = 1e-6
#     dt_check = dt/2.
#     while abs((distck3 - (r1 +r2))) > tol:
#         precollision = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_check)
#         distck3 = h3dist(rot2hyp(precollision[0:3]),rot2hyp(precollision[3:6]))

#         if (distck3 - (r1 +r2)) < 0:
#             dt_check = dt_check - tol
#         elif (distck3 - (r1 +r2)) > 0:
#             dt_check = dt_check + tol

#     postcollision = h3collision(precollision,params)

#     gs3simdatalist[step] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_check)
# else:
#     print("Not collided")

# startvecgs3sim = gs3simdatalist[step]

step += 1

while (step <= int(t_max/dt)):

    if distck3 < (r1 + r2):
        print("Collided")
        dt_check = dt/2.
        # step -= 1
        # nump = 0
        print("")
        print((distck3 - (r1 + r2)))
        if abs((distck3 - (r1 +r2))) > 1e-6:
            while abs((distck3 - (r1 +r2))) > 1e-6:# and nump < 100: # and dt-gs3cofact >= 1e-6 :
                precollision = gausss3(startvec=gs3simdatalist[step-2],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_check)
                distck3 = h3dist(rot2hyp(precollision[0:3]),rot2hyp(precollision[3:6]))
                # There are options on how to handle the time step scaling (addition or multiplication)
                # gs3cofact = gs3cofact * 1.01 #+ 1e-5
                # print(distck3)
                # print(dt-gs3cofact)
                # print((distck3 - (r1 +r2)))
                if (distck3 - (r1 +r2)) < 0:
                    dt_check = dt_check - 1e-6
                elif (distck3 - (r1 +r2)) > 0:
                    dt_check = dt_check + 1e-6
                print("dist")
                print(distck3)
                print("dist - r1 -r2")
                print((distck3 - (r1 +r2)))
                print("dt_check")
                print(dt_check)
                print("")
                # nump += 1
            postcollision = h3collision(precollision,params)
            # print(precollision)
            # print(postcollision)
            gs3simdatalist[step-1] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_check)
        else:
            postcollision = h3collision(gs3simdatalist[step-2],params)
            # print(precollision)
            # print(postcollision)
            gs3simdatalist[step-1] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
        step -= 1
    else:
        print("Not collided")
        gs3simdatalist[step] = gausss3(startvec=gs3simdatalist[step-1],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)


    startvecgs3sim = gs3simdatalist[step]

    # print("should be same as above")
    distck3 = h3dist(rot2hyp(startvecgs3sim[0:3]),rot2hyp(startvecgs3sim[3:6]))
    print("step {}".format(step))
    print(distck3)

    if step%int(1/dt)==0:
            print(step)
    step += 1


# Sim in H3
# np.save("h3_r_gausss1_sim_tmax10_dt1",gs1simdatalist)
# np.save("h3_r_gausss2_sim_tmax10_dt1",gs2simdatalist)
np.save("h3_r_gausss3_sim_tmax10_dt1",gs3simdatalist)


# Sim in S3
# np.save("s3_r_gausss1_sim_tmax10_dt0001",gs1simdatalist)
# np.save("s3_r_gausss2_sim_tmax10_dt0001",gs2simdatalist)
# np.save("s3_r_gausss3_sim_tmax10_dt0001",gs3simdatalist)

# data1 = np.load("h3_r_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("h3_r_gausss2_sim_tmax10_dt1.npy")
data3 = np.load("h3_r_gausss3_sim_tmax10_dt1.npy")

# data1 = np.load("s3_r_gausss1_sim_tmax10_dt0001.npy")
# data2 = np.load("s3_r_gausss2_sim_tmax10_dt0001.npy")
# data3 = np.load("s3_r_gausss3_sim_tmax10_dt0001.npy")




fig,ax=plt.subplots(1,1)

# distdata1 = np.zeros(np.shape(data1)[0])
# distdata2 = np.zeros(np.shape(data2)[0])
distdata3 = np.zeros(np.shape(data3)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
#     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
#     distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
#     counter += 1

counter = 0
for a in range(np.shape(data3)[0]):
    # distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
    # distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
    distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
    counter += 1

# ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")
# ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
# ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
ax.plot(t_arr,np.full(len(t_arr),1.))
ax.legend()
ax.set_title('Simulation Data')
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







