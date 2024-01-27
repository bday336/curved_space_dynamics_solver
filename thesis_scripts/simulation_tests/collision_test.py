# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from collision_function_bank import h3collision, s3collision, h3kvecproj
from test_system_bank import dynfunc_h3sim2ballcol, dynjac_h3sim2ballcol, dynfunc_s3sim2ballcol, dynjac_s3sim2ballcol

# Solver Setup

# Time array based on time step
dt = .001    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Time array based on number of steps
# nump = 10000    # Number of steps
# t_max = 10      # Total simulation time
# t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

# Simulation data container

# Sim H3
gs1simdatalist = np.zeros((t_arr.shape[0],12))
gs2simdatalist = np.zeros((t_arr.shape[0],12))
gs3simdatalist = np.zeros((t_arr.shape[0],12))

# Sim S3
# gs1simdatalist = np.zeros((t_arr.shape[0],12))
# gs2simdatalist = np.zeros((t_arr.shape[0],12))
# gs3simdatalist = np.zeros((t_arr.shape[0],12))

# Initial Data
v = .5      # Initial Velocity
x = 1.      # Spring Rest Length H3/S3
m1 = 1.      # Mass of point mass 1
m2 = 1.      # Mass of point mass 2
r1 = .5      # Radius of mass 1
r2 = .5      # Radisu of mass 2
params = [m1,m2,r1,r2]

# Collision tolerance
tol = 1e-8


# Sim col 1D H3
startvec = np.array([
    [2.,np.pi/2.,0.*np.pi/2.],[2.,np.pi/2.,2.*np.pi/2.],
    killingvech3([2.,np.pi/2.,0.*np.pi/2.],-v,"x"), killingvech3([2.,np.pi/2.,2.*np.pi/2.],v,"x")]).flatten()

# Sim col 2D H3
# startvec = np.array([
#     [2.,np.pi/2.,2.*np.pi/4.],[2.,np.pi/2.,0.*np.pi/4.],
#     killingvech3([2.,np.pi/2.,2.*np.pi/4.],-v,"y"), killingvech3([2.,np.pi/2.,0.*np.pi/4.],-v,"x")]).flatten()

# Sim col 1D S3
# startvec = np.array([
#     [np.pi/2.,np.pi/2.,np.pi/2.],[np.pi/2.,np.pi/2.,-np.pi/2.],
#     killingvecs3([np.pi/2.,np.pi/2.,np.pi/2.],-v,"vz"), killingvecs3([np.pi/2.,np.pi/2.,-np.pi/2.],v,"vz")]).flatten()

# Sim col 2D S3
# startvec = np.array([
#     [2.*np.pi/3. - .2,np.pi/2.,0.+.1],[2.*np.pi/3.,np.pi/2.,np.pi/2.-.1],
#     killingvecs3([2.*np.pi/3. - .2,np.pi/2.,0.+.1],v,"by"), killingvecs3([2.*np.pi/3.,np.pi/2.,np.pi/2.-.1],-v,"bx")]).flatten()



# Initial Data for step = 0

# Sim in H3
gs1simdatalist[0] = startvec.copy()
gs2simdatalist[0] = startvec.copy()
gs3simdatalist[0] = startvec.copy()

distck1 = h3dist(rot2hyp(gs1simdatalist[0][0:3]),rot2hyp(gs1simdatalist[0][3:6]))
distck2 = h3dist(rot2hyp(gs2simdatalist[0][0:3]),rot2hyp(gs2simdatalist[0][3:6]))
distck3 = h3dist(rot2hyp(gs3simdatalist[0][0:3]),rot2hyp(gs3simdatalist[0][3:6]))



# First Step
step1 = 1
step2 = 1
step3 = 1

# Sim in H3
gs1simdatalist[step1] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
gs2simdatalist[step2] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
gs3simdatalist[step3] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt) 

distck1 = h3dist(rot2hyp(gs1simdatalist[step1][0:3]),rot2hyp(gs1simdatalist[step1][3:6]))
distck2 = h3dist(rot2hyp(gs2simdatalist[step2][0:3]),rot2hyp(gs2simdatalist[step2][3:6]))
distck3 = h3dist(rot2hyp(gs3simdatalist[step3][0:3]),rot2hyp(gs3simdatalist[step3][3:6]))


step1 += 1
step2 += 1
step3 += 1

while (step1 <= int(t_max/dt) and step2<= int(t_max/dt) and step3 <= int(t_max/dt)):

    # H3 collision check

    if abs((distck3 - (r1 +r2))) < tol:
        print("Collided s3 at step {}".format(step3))
        dt_check = dt/2.
        # tol = 1e-6
        nump = 0

        if abs((distck3 - (r1 +r2))) > tol:
            print("Needed to find the collision position")
        
            dt_change = dt_check
            while abs((distck3 - (r1 +r2))) > tol and nump < 100:
                precollision = gausss3(startvec=gs3simdatalist[step3-2],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_change)
                distck3 = h3dist(rot2hyp(precollision[0:3]),rot2hyp(precollision[3:6]))

                # Still overlaping
                if (distck3 - (r1 +r2)) < 0:
                    dt_change = dt_change - dt_check #- dt_check*abs((distck3 - (r1 +r2))/(r1 +r2))/dt
                    print("dt_check minus")
                    print(dt_check)
                # Not overlapping
                elif (distck3 - (r1 +r2)) > 0:
                    dt_check = dt_check/2
                    dt_change = dt_change + dt_check #+ dt_check*abs((distck3 - (r1 +r2))/(r1 +r2))/dt
                    print("dt_check plus")
                    print(dt_check)
                print('separation')
                print(distck3 - (r1 +r2))
            

                nump += 1
            prepart1x = h3kvecproj(precollision[0:3],precollision[6:9],"x")
            prepart1y = h3kvecproj(precollision[0:3],precollision[6:9],"y")
            prepart1z = h3kvecproj(precollision[0:3],precollision[6:9],"z")

            prepart2x = h3kvecproj(precollision[3:6],precollision[9:12],"x")
            prepart2y = h3kvecproj(precollision[3:6],precollision[9:12],"y")
            prepart2z = h3kvecproj(precollision[3:6],precollision[9:12],"z")

            print("Pre Projection of particle1")
            print([prepart1x,prepart1y,prepart1z])
            print("Pre Projection of particle2")
            print([prepart2x,prepart2y,prepart2z])
            print("")
            print("Total Pre Projection")
            print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
            print("")

            postcollision = h3collision(precollision,params)

            postpart1x = h3kvecproj(postcollision[0:3],postcollision[6:9],"x")
            postpart1y = h3kvecproj(postcollision[0:3],postcollision[6:9],"y")
            postpart1z = h3kvecproj(postcollision[0:3],postcollision[6:9],"z")

            postpart2x = h3kvecproj(postcollision[3:6],postcollision[9:12],"x")
            postpart2y = h3kvecproj(postcollision[3:6],postcollision[9:12],"y")
            postpart2z = h3kvecproj(postcollision[3:6],postcollision[9:12],"z")

            print("Post Projection of particle1")
            print([postpart1x,postpart1y,postpart1z])
            print("Post Projection of particle2")
            print([postpart2x,postpart2y,postpart2z])
            print("")
            print("Total Post Projection")
            print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z]))
            print("")
            print("Difference between pre and post")
            print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z])-np.array([prepart1x,prepart1y,prepart1z])-np.array([prepart2x,prepart2y,prepart2z]))
            print("")

            gs3simdatalist[step3-1] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_change)
            step3 -= 1
            break
        # If the collision happens to be within tolerance
        else:
            print("Collision occurred within tolerance")

            prepart1x = h3kvecproj(gs3simdatalist[step3-1][0:3],gs3simdatalist[step3-1][6:9],"x")
            prepart1y = h3kvecproj(gs3simdatalist[step3-1][0:3],gs3simdatalist[step3-1][6:9],"y")
            prepart1z = h3kvecproj(gs3simdatalist[step3-1][0:3],gs3simdatalist[step3-1][6:9],"z")

            prepart2x = h3kvecproj(gs3simdatalist[step3-1][3:6],gs3simdatalist[step3-1][9:12],"x")
            prepart2y = h3kvecproj(gs3simdatalist[step3-1][3:6],gs3simdatalist[step3-1][9:12],"y")
            prepart2z = h3kvecproj(gs3simdatalist[step3-1][3:6],gs3simdatalist[step3-1][9:12],"z")

            print("Pre Projection of particle1")
            print([prepart1x,prepart1y,prepart1z])
            print("Pre Projection of particle2")
            print([prepart2x,prepart2y,prepart2z])
            print("")
            print("Total Pre Projection")
            print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
            print("")

            postcollision = h3collision(gs3simdatalist[step3-1],params)

            postpart1x = h3kvecproj(postcollision[0:3],postcollision[6:9],"x")
            postpart1y = h3kvecproj(postcollision[0:3],postcollision[6:9],"y")
            postpart1z = h3kvecproj(postcollision[0:3],postcollision[6:9],"z")

            postpart2x = h3kvecproj(postcollision[3:6],postcollision[9:12],"x")
            postpart2y = h3kvecproj(postcollision[3:6],postcollision[9:12],"y")
            postpart2z = h3kvecproj(postcollision[3:6],postcollision[9:12],"z")

            print("Post Projection of particle1")
            print([postpart1x,postpart1y,postpart1z])
            print("Post Projection of particle2")
            print([postpart2x,postpart2y,postpart2z])
            print("")
            print("Total Post Projection")
            print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z]))
            print("")
            print("Difference between pre and post")
            print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z])-np.array([prepart1x,prepart1y,prepart1z])-np.array([prepart2x,prepart2y,prepart2z]))
            print("")

            gs3simdatalist[step3] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
    else:
        # print("Not collided")
        gs3simdatalist[step3] = gausss3(startvec=gs3simdatalist[step3-1],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)


    # Sim in H3
    startvecgs1sim = gs1simdatalist[step1]
    startvecgs2sim = gs2simdatalist[step2]
    startvecgs3sim = gs3simdatalist[step3]

    distck1 = h3dist(rot2hyp(startvecgs1sim[0:3]),rot2hyp(startvecgs1sim[3:6]))
    distck2 = h3dist(rot2hyp(startvecgs2sim[0:3]),rot2hyp(startvecgs2sim[3:6]))
    distck3 = h3dist(rot2hyp(startvecgs3sim[0:3]),rot2hyp(startvecgs3sim[3:6]))


    # print("dist3")
    # print(distck3)
    # print("dist3 - r1 -r2")
    # print((distck3 - (r1 +r2)))


    # print("step")
    # print(step3)
    # print("dist3")
    # print(distck3)
    # print("dist3 - r1 -r2")
    # print((distck3 - (r1 +r2)))
    # if step1%int(1/dt)==0:
    #         print(step1)
    # if step2%int(1/dt)==0:
    #         print(step2)
    # if step3%int(1/dt)==0:
    #         print(step3)
    step1 += 1
    step2 += 1
    step3 += 1


# # Sim in H3
# np.save("h3_1d_col_gausss1_sim_tmax10_dt1",gs1simdatalist)
# np.save("h3_1d_col_gausss2_sim_tmax10_dt1",gs2simdatalist)
np.save("h3_1d_col_gausss3_sim_tmax10_dt1",gs3simdatalist)

# data1 = np.load("h3_1d_col_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("h3_1d_col_gausss2_sim_tmax10_dt1.npy")
data3 = np.load("h3_1d_col_gausss3_sim_tmax10_dt1.npy")

#Sim in S3
# np.save("s3_1d_col_gausss1_sim_tmax10_dt1",gs1simdatalist)
# np.save("s3_1d_col_gausss2_sim_tmax10_dt1",gs2simdatalist)
# np.save("s3_1d_col_gausss3_sim_tmax10_dt1",gs3simdatalist)

# data1 = np.load("s3_1d_col_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("s3_1d_col_gausss2_sim_tmax10_dt1.npy")
# data3 = np.load("s3_1d_col_gausss3_sim_tmax10_dt1.npy")




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
ax.plot(t_arr,distdata3,'b',marker='o',label = "Gauss s3")
ax.plot(t_arr,np.full(len(t_arr),1.))
ax.legend()
ax.set_title('Simulation Data')
# ax.set_xlim(2.95,3.05)
# ax.set_ylim(.95,1.05)
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







