# Generic Simulation Script

import numpy as np
import matplotlib.pyplot as plt
from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3,boostxh3
from integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from collision_function_bank import h3collision, s3collision, h3kvecproj, s3kvecproj,transform2origh3col
from test_system_bank import dynfunc_h3sim2ballcol, dynjac_h3sim2ballcol, dynfunc_s3sim2ballcol, dynjac_s3sim2ballcol

# Solver Setup

# Time array based on time step
dt = .01    # Number of steps
t_max = 10      # Total simulation time
t_arr = np.arange(0.,t_max+dt,dt)

# Time array based on number of steps
# nump = 10000    # Number of steps
# t_max = 10      # Total simulation time
# t_arr, dt= np.linspace(0.,t_max,nump,retstep=True)

# Simulation data container

# Sim H3
# gs1simdatalist = np.zeros((t_arr.shape[0],12))
# gs2simdatalist = np.zeros((t_arr.shape[0],12))
gs3simdatalist = np.zeros((t_arr.shape[0],12))

# Sim S3
# gs1simdatalist = np.zeros((t_arr.shape[0],12))
# gs2simdatalist = np.zeros((t_arr.shape[0],12))
# gs3simdatalist = np.zeros((t_arr.shape[0],12))

# Initial Data
v = .3      # Initial Velocity
x = 1.      # Spring Rest Length H3/S3
m1 = 1.      # Mass of point mass 1
m2 = 1.      # Mass of point mass 2
# r1 = .5      # Radius of mass 1 H3
# r2 = .5      # Radisu of mass 2 H3
r1 = .5      # Radius of mass 1 S3
r2 = .5     # Radisu of mass 2 S3
params = [m1,m2,r1,r2]

# Collision tolerance
tol = 1e-20


# Sim col 1D H3
# startvec = np.array([
#     [2.,np.pi/2.,0.*np.pi/2.],[2.,np.pi/2.,2.*np.pi/2.],
#     killingvech3([2.,np.pi/2.,0.*np.pi/2.],-v,"x"), killingvech3([2.,np.pi/2.,2.*np.pi/2.],v,"x")]).flatten()

# Sim col 2D H3
# startvec = np.array([
#     [2.,np.pi/2.,2.*np.pi/4.],[2.,np.pi/2.,0.*np.pi/4.],
#     killingvech3([2.,np.pi/2.,2.*np.pi/4.],-v,"y"), killingvech3([2.,np.pi/2.,0.*np.pi/4.],-v,"x")]).flatten()
# startvec = np.array([
#     [2.,np.pi/2.,2.*np.pi/4.],[2.,np.pi/2.,0.*np.pi/4.],
#     [-v,0,0],[-v,0,0]]).flatten()
## Angled version
### V
# startvec = np.array([
#     [2.,np.pi/2.,7.*np.pi/16.+np.pi/2.],[2.,np.pi/2.,25.*np.pi/16.+np.pi/2.],
#     [-.5,0,0],[-.5,0,0]]).flatten()
## < Used this one to test
# startvec = np.array([
#     [2.,np.pi/2.,7.*np.pi/16.],[2.,np.pi/2.,25.*np.pi/16.],
#     [-.5,0,0],[-.5,0,0]]).flatten()

# Sim col 1D S3
# startvec = np.array([
#     [np.pi/2.,np.pi/2.,np.pi/2.],[np.pi/2.,np.pi/2.,-np.pi/2.],
#     killingvecs3([np.pi/2.,np.pi/2.,np.pi/2.],-v,"vz"), killingvecs3([np.pi/2.,np.pi/2.,-np.pi/2.],v,"vz")]).flatten()

# Sim col 2D S3
# startvec = np.array([
#     [2.*np.pi/3. - .2,np.pi/2.,0.+.1],[2.*np.pi/3.,np.pi/2.,np.pi/2.-.1],
#     killingvecs3([2.*np.pi/3. - .2,np.pi/2.,0.+.1],v,"by"), killingvecs3([2.*np.pi/3.,np.pi/2.,np.pi/2.-.1],-v,"bx")]).flatten()
## V
startvec = np.array([
    [np.pi/3.,np.pi/2.,-np.pi/2.],[np.pi/3.,np.pi/2.,np.pi/2.],
    [-.1,0,.2],[-.1,0,-.2]]).flatten()



# Sim in H3
# # gs1simdatalist[0] = startvec.copy()
# # gs2simdatalist[0] = startvec.copy()
# gs3simdatalist[0] = startvec.copy()

# # distck1 = h3dist(rot2hyp(gs1simdatalist[0][0:3]),rot2hyp(gs1simdatalist[0][3:6]))
# # distck2 = h3dist(rot2hyp(gs2simdatalist[0][0:3]),rot2hyp(gs2simdatalist[0][3:6]))
# distck3 = h3dist(rot2hyp(gs3simdatalist[0][0:3]),rot2hyp(gs3simdatalist[0][3:6]))


# Sim in S3
# gs1simdatalist[0] = startvec.copy()
# gs2simdatalist[0] = startvec.copy()
gs3simdatalist[0] = startvec.copy()

# distck1 = r4dist(rot2r4(gs1simdatalist[0][0:3]),rot2r4(gs1simdatalist[0][3:6]))
# distck2 = r4dist(rot2r4(gs2simdatalist[0][0:3]),rot2r4(gs2simdatalist[0][3:6]))
distck3 = r4dist(rot2r4(gs3simdatalist[0][0:3]),rot2r4(gs3simdatalist[0][3:6]))

# First Step
# step1 = 0
# step2 = 0
step3 = 0

# Sim in H3
# # gs1simdatalist[step1] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
# # gs2simdatalist[step2] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
# gs3simdatalist[step3] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt) 

# # distck1 = h3dist(rot2hyp(gs1simdatalist[step1][0:3]),rot2hyp(gs1simdatalist[step1][3:6]))
# # distck2 = h3dist(rot2hyp(gs2simdatalist[step2][0:3]),rot2hyp(gs2simdatalist[step2][3:6]))
# distck3 = h3dist(rot2hyp(gs3simdatalist[step3][0:3]),rot2hyp(gs3simdatalist[step3][3:6]))

# Sim in S3
# gs1simdatalist[step1] = gausss1(startvec=startvec,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)
# gs2simdatalist[step2] = gausss2(startvec=startvec,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)
gs3simdatalist[step3] = gausss3(startvec=startvec,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt) 

# distck1 = r4dist(rot2r4(gs1simdatalist[step1][0:3]),rot2r4(gs1simdatalist[step1][3:6]))
# distck2 = r4dist(rot2r4(gs2simdatalist[step2][0:3]),rot2r4(gs2simdatalist[step2][3:6]))
distck3 = r4dist(rot2r4(gs3simdatalist[step3][0:3]),rot2r4(gs3simdatalist[step3][3:6]))

#### We assume that the system does not start in a collided state

# step1 += 1
# step2 += 1
step3 += 1

while (step3 <= int(t_max/dt)):
#  while (step1 < int(t_max/dt) and step2< int(t_max/dt) and step3 < int(t_max/dt)):

    # H3 collision check

    # if abs((distck1 - (r1 +r2))) < tol or (distck1 - (r1 +r2)) < 0 :
    #     print("Collided s1")
    #     dt_check = dt/2.
    #     # tol = 1e-6
    #     nump = 0

    #     if abs((distck1 - (r1 +r2))) > tol:
    #         print("Needed to find the collision position")

    #         dt_change = dt_check
    #         while abs((distck1 - (r1 +r2))) > tol and nump < 100:
    #             precollision = gausss1(startvec=gs1simdatalist[step1-2],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_change)
    #             distck1 = h3dist(rot2hyp(precollision[0:3]),rot2hyp(precollision[3:6]))
    #             if (distck1 - (r1 +r2)) < 0:
    #                 dt_change = dt_change - dt_check
    #                 print("dt_check minus")
    #                 print(dt_check)
    #             elif (distck1 - (r1 +r2)) > 0:
    #                 dt_check = dt_check/2
    #                 dt_change = dt_change + dt_check
    #                 print("dt_check plus")
    #                 print(dt_check)

    #         prepart1bx, prepart1bxvec = h3kvecproj(precollision[0:3],precollision[6:9],"bx")
    #         prepart1by, prepart1byvec = h3kvecproj(precollision[0:3],precollision[6:9],"by")
    #         prepart1bz, prepart1bzvec = h3kvecproj(precollision[0:3],precollision[6:9],"bz")
    #         prepart1vx, prepart1vxvec = h3kvecproj(precollision[0:3],precollision[6:9],"wx")
    #         prepart1vy, prepart1vyvec = h3kvecproj(precollision[0:3],precollision[6:9],"wy")
    #         prepart1vz, prepart1vzvec = h3kvecproj(precollision[0:3],precollision[6:9],"wz")

    #         prepart2bx, prepart2bxvec = h3kvecproj(precollision[3:6],precollision[9:12],"bx")
    #         prepart2by, prepart2byvec = h3kvecproj(precollision[3:6],precollision[9:12],"by")
    #         prepart2bz, prepart2bzvec = h3kvecproj(precollision[3:6],precollision[9:12],"bz")
    #         prepart2vx, prepart2vxvec = h3kvecproj(precollision[3:6],precollision[9:12],"wx")
    #         prepart2vy, prepart2vyvec = h3kvecproj(precollision[3:6],precollision[9:12],"wy")
    #         prepart2vz, prepart2vzvec = h3kvecproj(precollision[3:6],precollision[9:12],"wz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])+np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         # print("")

    #         postcollision = h3collision(precollision,params)

    #         postpart1bx, postpart1bxvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by, postpart1byvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz, postpart1bzvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx, postpart1vxvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wx")
    #         postpart1vy, postpart1vyvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wy")
    #         postpart1vz, postpart1vzvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wz")

    #         postpart2bx, postpart2bxvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by, postpart2byvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz, postpart2bzvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx, postpart2vxvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wx")
    #         postpart2vy, postpart2vyvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wy")
    #         postpart2vz, postpart2vzvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])+np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs1simdatalist[step1-1] = gausss1(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_change)
    #         step1 -= 1
    #     # If the collision happens to be within tolerance
    #     else:
    #         print("Collision occurred within tolerance")

    #         prepart1bx = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"bx")
    #         prepart1by = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"by")
    #         prepart1bz = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"bz")
    #         prepart1vx = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"wx")
    #         prepart1vy = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"wy")
    #         prepart1vz = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"wz")

    #         prepart2bx = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"bx")
    #         prepart2by = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"by")
    #         prepart2bz = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"bz")
    #         prepart2vx = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"wx")
    #         prepart2vy = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"wy")
    #         prepart2vz = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"wz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
    #         # print("")

    #         postcollision = h3collision(gs1simdatalist[step1-1],params)

    #         postpart1bx = h3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by = h3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz = h3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx = h3kvecproj(postcollision[0:3],postcollision[6:9],"wx")
    #         postpart1vy = h3kvecproj(postcollision[0:3],postcollision[6:9],"wy")
    #         postpart1vz = h3kvecproj(postcollision[0:3],postcollision[6:9],"wz")

    #         postpart2bx = h3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by = h3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz = h3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx = h3kvecproj(postcollision[3:6],postcollision[9:12],"wx")
    #         postpart2vy = h3kvecproj(postcollision[3:6],postcollision[9:12],"wy")
    #         postpart2vz = h3kvecproj(postcollision[3:6],postcollision[9:12],"wz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs1simdatalist[step1] = gausss1(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
    # else:
    #     # print("Not collided")
    #     gs1simdatalist[step1] = gausss3(startvec=gs1simdatalist[step1-1],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)



    # if abs((distck2 - (r1 +r2))) < tol or (distck2 - (r1 +r2)) < 0:
    #     print("Collided s2")
    #     dt_check = dt/2.
    #     # tol = 1e-6
    #     nump = 0

    #     if abs((distck2 - (r1 +r2))) > tol:
    #         print("Needed to find the collision position")

    #         dt_change = dt_check
    #         while abs((distck2 - (r1 +r2))) > tol and nump < 100:
    #             precollision = gausss2(startvec=gs2simdatalist[step2-2],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_change)
    #             distck2 = h3dist(rot2hyp(precollision[0:3]),rot2hyp(precollision[3:6]))
    #             if (distck2 - (r1 +r2)) < 0:
    #                 dt_change = dt_change - dt_check
    #                 # print("dt_check minus")
    #                 # print(dt_check)
    #             elif (distck2 - (r1 +r2)) > 0:
    #                 dt_check = dt_check/2
    #                 dt_change = dt_change + dt_check
    #                 # print("dt_check plus")
    #                 # print(dt_check)
    #             nump += 1
                
    #         prepart1bx = h3kvecproj(precollision[0:3],precollision[6:9],"bx")
    #         prepart1by = h3kvecproj(precollision[0:3],precollision[6:9],"by")
    #         prepart1bz = h3kvecproj(precollision[0:3],precollision[6:9],"bz")
    #         prepart1vx = h3kvecproj(precollision[0:3],precollision[6:9],"wx")
    #         prepart1vy = h3kvecproj(precollision[0:3],precollision[6:9],"wy")
    #         prepart1vz = h3kvecproj(precollision[0:3],precollision[6:9],"wz")

    #         prepart2bx = h3kvecproj(precollision[3:6],precollision[9:12],"bx")
    #         prepart2by = h3kvecproj(precollision[3:6],precollision[9:12],"by")
    #         prepart2bz = h3kvecproj(precollision[3:6],precollision[9:12],"bz")
    #         prepart2vx = h3kvecproj(precollision[3:6],precollision[9:12],"wx")
    #         prepart2vy = h3kvecproj(precollision[3:6],precollision[9:12],"wy")
    #         prepart2vz = h3kvecproj(precollision[3:6],precollision[9:12],"wz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])+np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         # print("")

    #         postcollision = h3collision(precollision,params)

    #         postpart1bx = h3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by = h3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz = h3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx = h3kvecproj(postcollision[0:3],postcollision[6:9],"wx")
    #         postpart1vy = h3kvecproj(postcollision[0:3],postcollision[6:9],"wy")
    #         postpart1vz = h3kvecproj(postcollision[0:3],postcollision[6:9],"wz")

    #         postpart2bx = h3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by = h3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz = h3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx = h3kvecproj(postcollision[3:6],postcollision[9:12],"wx")
    #         postpart2vy = h3kvecproj(postcollision[3:6],postcollision[9:12],"wy")
    #         postpart2vz = h3kvecproj(postcollision[3:6],postcollision[9:12],"wz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])+np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs2simdatalist[step2-1] = gausss2(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_change)
    #         step2 -= 1
    #     # If the collision happens to be within tolerance
    #     else:
    #         print("Collision occurred within tolerance")

    #         prepart1bx = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"bx")
    #         prepart1by = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"by")
    #         prepart1bz = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"bz")
    #         prepart1vx = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"wx")
    #         prepart1vy = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"wy")
    #         prepart1vz = h3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"wz")

    #         prepart2bx = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"bx")
    #         prepart2by = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"by")
    #         prepart2bz = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"bz")
    #         prepart2vx = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"wx")
    #         prepart2vy = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"wy")
    #         prepart2vz = h3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"wz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
    #         # print("")

    #         postcollision = h3collision(gs1simdatalist[step1-1],params)

    #         postpart1bx = h3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by = h3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz = h3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx = h3kvecproj(postcollision[0:3],postcollision[6:9],"wx")
    #         postpart1vy = h3kvecproj(postcollision[0:3],postcollision[6:9],"wy")
    #         postpart1vz = h3kvecproj(postcollision[0:3],postcollision[6:9],"wz")

    #         postpart2bx = h3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by = h3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz = h3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx = h3kvecproj(postcollision[3:6],postcollision[9:12],"wx")
    #         postpart2vy = h3kvecproj(postcollision[3:6],postcollision[9:12],"wy")
    #         postpart2vz = h3kvecproj(postcollision[3:6],postcollision[9:12],"wz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs2simdatalist[step2] = gausss2(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
    # else:
    #     # print("Not collided")
    #     gs2simdatalist[step2] = gausss2(startvec=gs2simdatalist[step2-1],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)



    # if abs((distck3 - (r1 +r2))) < tol or (distck3 - (r1 +r2)) < 0:
    #     print("Collided s3 at step {}".format(step3))
    #     dt_check = dt/2.
    #     # tol = 1e-6
    #     nump = 0

    #     if abs((distck3 - (r1 +r2))) > tol:
    #         print("Needed to find the collision position")
        
    #         dt_change = dt_check
    #         while abs((distck3 - (r1 +r2))) > tol and nump < 100:
    #             precollision = gausss3(startvec=gs3simdatalist[step3-2],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt_change)
    #             distck3 = h3dist(rot2hyp(precollision[0:3]),rot2hyp(precollision[3:6]))

    #             # Still overlaping
    #             if (distck3 - (r1 +r2)) < 0:
    #                 dt_change = dt_change - dt_check
    #                 # print("dt_check minus")
    #                 # print(dt_check)
    #             # Not overlapping
    #             elif (distck3 - (r1 +r2)) > 0:
    #                 dt_check = dt_check/2
    #                 dt_change = dt_change + dt_check
    #                 # print("dt_check plus")
    #                 # print(dt_check)
    #             # print('separation')
    #             # print(distck3 - (r1 +r2))

    #             nump += 1

    #         prepart1bx, prepart1bxvec = h3kvecproj(precollision[0:3],precollision[6:9],"bx")
    #         prepart1by, prepart1byvec = h3kvecproj(precollision[0:3],precollision[6:9],"by")
    #         prepart1bz, prepart1bzvec = h3kvecproj(precollision[0:3],precollision[6:9],"bz")
    #         prepart1vx, prepart1vxvec = h3kvecproj(precollision[0:3],precollision[6:9],"wx")
    #         prepart1vy, prepart1vyvec = h3kvecproj(precollision[0:3],precollision[6:9],"wy")
    #         prepart1vz, prepart1vzvec = h3kvecproj(precollision[0:3],precollision[6:9],"wz")

    #         prepart2bx, prepart2bxvec = h3kvecproj(precollision[3:6],precollision[9:12],"bx")
    #         prepart2by, prepart2byvec = h3kvecproj(precollision[3:6],precollision[9:12],"by")
    #         prepart2bz, prepart2bzvec = h3kvecproj(precollision[3:6],precollision[9:12],"bz")
    #         prepart2vx, prepart2vxvec = h3kvecproj(precollision[3:6],precollision[9:12],"wx")
    #         prepart2vy, prepart2vyvec = h3kvecproj(precollision[3:6],precollision[9:12],"wy")
    #         prepart2vz, prepart2vzvec = h3kvecproj(precollision[3:6],precollision[9:12],"wz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])+np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         # print("")

    #         postcollision = h3collision(precollision,params)

    #         postpart1bx, postpart1bxvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by, postpart1byvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz, postpart1bzvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx, postpart1vxvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wx")
    #         postpart1vy, postpart1vyvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wy")
    #         postpart1vz, postpart1vzvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wz")

    #         postpart2bx, postpart2bxvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by, postpart2byvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz, postpart2bzvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx, postpart2vxvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wx")
    #         postpart2vy, postpart2vyvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wy")
    #         postpart2vz, postpart2vzvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])+np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs3simdatalist[step3-1] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_change)
    #         step3 -= 1
    #         # break
    #     # If the collision happens to be within tolerance
    #     else:
    #         print("Collision occurred within tolerance")

    #         prepart1bx, prepart1bxvec = h3kvecproj(gs3simdatalist[step1-1][0:3],gs3simdatalist[step1-1][6:9],"bx")
    #         prepart1by, prepart1byvec = h3kvecproj(gs3simdatalist[step1-1][0:3],gs3simdatalist[step1-1][6:9],"by")
    #         prepart1bz, prepart1bzvec = h3kvecproj(gs3simdatalist[step1-1][0:3],gs3simdatalist[step1-1][6:9],"bz")
    #         prepart1vx, prepart1vxvec = h3kvecproj(gs3simdatalist[step1-1][0:3],gs3simdatalist[step1-1][6:9],"wx")
    #         prepart1vy, prepart1vyvec = h3kvecproj(gs3simdatalist[step1-1][0:3],gs3simdatalist[step1-1][6:9],"wy")
    #         prepart1vz, prepart1vzvec = h3kvecproj(gs3simdatalist[step1-1][0:3],gs3simdatalist[step1-1][6:9],"wz")

    #         prepart2bx, prepart2bxvec = h3kvecproj(gs3simdatalist[step1-1][3:6],gs3simdatalist[step1-1][9:12],"bx")
    #         prepart2by, prepart2byvec = h3kvecproj(gs3simdatalist[step1-1][3:6],gs3simdatalist[step1-1][9:12],"by")
    #         prepart2bz, prepart2bzvec = h3kvecproj(gs3simdatalist[step1-1][3:6],gs3simdatalist[step1-1][9:12],"bz")
    #         prepart2vx, prepart2vxvec = h3kvecproj(gs3simdatalist[step1-1][3:6],gs3simdatalist[step1-1][9:12],"wx")
    #         prepart2vy, prepart2vyvec = h3kvecproj(gs3simdatalist[step1-1][3:6],gs3simdatalist[step1-1][9:12],"wy")
    #         prepart2vz, prepart2vzvec = h3kvecproj(gs3simdatalist[step1-1][3:6],gs3simdatalist[step1-1][9:12],"wz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
    #         # print("")

    #         postcollision = h3collision(gs1simdatalist[step1-1],params)

    #         postpart1bx, postpart1bxvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by, postpart1byvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz, postpart1bzvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx, postpart1vxvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wx")
    #         postpart1vy, postpart1vyvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wy")
    #         postpart1vz, postpart1vzvec = h3kvecproj(postcollision[0:3],postcollision[6:9],"wz")

    #         postpart2bx, postpart2bxvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by, postpart2byvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz, postpart2bzvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx, postpart2vxvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wx")
    #         postpart2vy, postpart2vyvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wy")
    #         postpart2vz, postpart2vzvec = h3kvecproj(postcollision[3:6],postcollision[9:12],"wz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs3simdatalist[step3] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)
    # else:
    #     # print("Not collided")
    #     gs3simdatalist[step3] = gausss3(startvec=gs3simdatalist[step3-1],params=params,dynfunc=dynfunc_h3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt)

    # # S3 collision check

    # if abs((distck1 - (r1 +r2))) < tol or (distck1 - (r1 +r2)) < 0:
    #     print("Collided s1")
    #     dt_check = dt/2.
    #     # tol = 1e-6
    #     nump = 0

    #     if abs((distck1 - (r1 +r2))) > tol:
    #         print("Needed to find the collision position")
        
    #         dt_change = dt_check
    #         while abs((distck1 - (r1 +r2))) > tol and nump < 100:
    #             precollision = gausss1(startvec=gs1simdatalist[step1-2],params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt_change)
    #             distck1 = r4dist(rot2r4(precollision[0:3]),rot2r4(precollision[3:6]))
    #             if (distck1 - (r1 +r2)) < 0:
    #                 dt_change = dt_change - dt_check
    #                 # print("dt_check minus")
    #                 # print(dt_check)
    #             elif (distck1 - (r1 +r2)) > 0:
    #                 dt_check = dt_check/2
    #                 dt_change = dt_change + dt_check
    #                 # print("dt_check plus")
    #                 # print(dt_check)
    #             nump += 1

    #         prepart1bx = s3kvecproj(precollision[0:3],precollision[6:9],"bx")
    #         prepart1by = s3kvecproj(precollision[0:3],precollision[6:9],"by")
    #         prepart1bz = s3kvecproj(precollision[0:3],precollision[6:9],"bz")
    #         prepart1vx = s3kvecproj(precollision[0:3],precollision[6:9],"vx")
    #         prepart1vy = s3kvecproj(precollision[0:3],precollision[6:9],"vy")
    #         prepart1vz = s3kvecproj(precollision[0:3],precollision[6:9],"vz")

    #         prepart2bx = s3kvecproj(precollision[3:6],precollision[9:12],"bx")
    #         prepart2by = s3kvecproj(precollision[3:6],precollision[9:12],"by")
    #         prepart2bz = s3kvecproj(precollision[3:6],precollision[9:12],"bz")
    #         prepart2vx = s3kvecproj(precollision[3:6],precollision[9:12],"vx")
    #         prepart2vy = s3kvecproj(precollision[3:6],precollision[9:12],"vy")
    #         prepart2vz = s3kvecproj(precollision[3:6],precollision[9:12],"vz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])+np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         # print("")

    #         postcollision = s3collision(precollision,params)

    #         postpart1bx = s3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by = s3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz = s3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx = s3kvecproj(postcollision[0:3],postcollision[6:9],"vx")
    #         postpart1vy = s3kvecproj(postcollision[0:3],postcollision[6:9],"vy")
    #         postpart1vz = s3kvecproj(postcollision[0:3],postcollision[6:9],"vz")

    #         postpart2bx = s3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by = s3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz = s3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx = s3kvecproj(postcollision[3:6],postcollision[9:12],"vx")
    #         postpart2vy = s3kvecproj(postcollision[3:6],postcollision[9:12],"vy")
    #         postpart2vz = s3kvecproj(postcollision[3:6],postcollision[9:12],"vz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])+np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs1simdatalist[step1-1] = gausss1(startvec=postcollision,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_change)
    #         step1 -= 1
    #         # break
    #     # If the collision happens to be within tolerance
    #     else:
    #         print("Collision occurred within tolerance")

    #         prepart1x = s3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"x")
    #         prepart1y = s3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"y")
    #         prepart1z = s3kvecproj(gs1simdatalist[step1-1][0:3],gs1simdatalist[step1-1][6:9],"z")

    #         prepart2x = s3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"x")
    #         prepart2y = s3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"y")
    #         prepart2z = s3kvecproj(gs1simdatalist[step1-1][3:6],gs1simdatalist[step1-1][9:12],"z")

    #         print("Pre Projection of particle1")
    #         print([prepart1x,prepart1y,prepart1z])
    #         print("Pre Projection of particle2")
    #         print([prepart2x,prepart2y,prepart2z])
    #         print("")
    #         print("Total Pre Projection")
    #         print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
    #         print("")

    #         postcollision = s3collision(gs1simdatalist[step1-1],params)

    #         postpart1x = s3kvecproj(postcollision[0:3],postcollision[6:9],"x")
    #         postpart1y = s3kvecproj(postcollision[0:3],postcollision[6:9],"y")
    #         postpart1z = s3kvecproj(postcollision[0:3],postcollision[6:9],"z")

    #         postpart2x = s3kvecproj(postcollision[3:6],postcollision[9:12],"x")
    #         postpart2y = s3kvecproj(postcollision[3:6],postcollision[9:12],"y")
    #         postpart2z = s3kvecproj(postcollision[3:6],postcollision[9:12],"z")

    #         print("Post Projection of particle1")
    #         print([postpart1x,postpart1y,postpart1z])
    #         print("Post Projection of particle2")
    #         print([postpart2x,postpart2y,postpart2z])
    #         print("")
    #         print("Total Post Projection")
    #         print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z]))
    #         print("")
    #         print("Difference between pre and post")
    #         print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z])-np.array([prepart1x,prepart1y,prepart1z])-np.array([prepart2x,prepart2y,prepart2z]))
    #         print("")

    #         gs1simdatalist[step1] = gausss1(startvec=postcollision,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)
    # else:
    #     # print("Not collided")
    #     gs1simdatalist[step1] = gausss1(startvec=gs1simdatalist[step1-1],params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)


    # if abs((distck2 - (r1 +r2))) < tol or (distck2 - (r1 +r2)) < 0:
    #     print("Collided s2")
    #     dt_check = dt/2.
    #     # tol = 1e-6
    #     nump = 0

    #     if abs((distck2 - (r1 +r2))) > tol:
    #         print("Needed to find the collision position")
        
    #         dt_change = dt_check
    #         while abs((distck2 - (r1 +r2))) > tol and nump < 100:
    #             precollision = gausss2(startvec=gs2simdatalist[step2-2],params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt_change)
    #             distck2 = r4dist(rot2r4(precollision[0:3]),rot2r4(precollision[3:6]))
    #             if (distck2 - (r1 +r2)) < 0:
    #                 dt_change = dt_change - dt_check
    #                 # print("dt_check minus")
    #                 # print(dt_check)
    #             elif (distck2 - (r1 +r2)) > 0:
    #                 dt_check = dt_check/2
    #                 dt_change = dt_change + dt_check
    #                 # print("dt_check plus")
    #                 # print(dt_check)
    #             nump += 1

    #         prepart1bx = s3kvecproj(precollision[0:3],precollision[6:9],"bx")
    #         prepart1by = s3kvecproj(precollision[0:3],precollision[6:9],"by")
    #         prepart1bz = s3kvecproj(precollision[0:3],precollision[6:9],"bz")
    #         prepart1vx = s3kvecproj(precollision[0:3],precollision[6:9],"vx")
    #         prepart1vy = s3kvecproj(precollision[0:3],precollision[6:9],"vy")
    #         prepart1vz = s3kvecproj(precollision[0:3],precollision[6:9],"vz")

    #         prepart2bx = s3kvecproj(precollision[3:6],precollision[9:12],"bx")
    #         prepart2by = s3kvecproj(precollision[3:6],precollision[9:12],"by")
    #         prepart2bz = s3kvecproj(precollision[3:6],precollision[9:12],"bz")
    #         prepart2vx = s3kvecproj(precollision[3:6],precollision[9:12],"vx")
    #         prepart2vy = s3kvecproj(precollision[3:6],precollision[9:12],"vy")
    #         prepart2vz = s3kvecproj(precollision[3:6],precollision[9:12],"vz")

    #         print("Pre Projection of particle1")
    #         print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
    #         print("Pre Projection of particle2")
    #         print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
    #         print("")
    #         # print("Total Pre Projection")
    #         # print(np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])+np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         # print("")

    #         postcollision = s3collision(precollision,params)

    #         postpart1bx = s3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
    #         postpart1by = s3kvecproj(postcollision[0:3],postcollision[6:9],"by")
    #         postpart1bz = s3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
    #         postpart1vx = s3kvecproj(postcollision[0:3],postcollision[6:9],"vx")
    #         postpart1vy = s3kvecproj(postcollision[0:3],postcollision[6:9],"vy")
    #         postpart1vz = s3kvecproj(postcollision[0:3],postcollision[6:9],"vz")

    #         postpart2bx = s3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
    #         postpart2by = s3kvecproj(postcollision[3:6],postcollision[9:12],"by")
    #         postpart2bz = s3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
    #         postpart2vx = s3kvecproj(postcollision[3:6],postcollision[9:12],"vx")
    #         postpart2vy = s3kvecproj(postcollision[3:6],postcollision[9:12],"vy")
    #         postpart2vz = s3kvecproj(postcollision[3:6],postcollision[9:12],"vz")

    #         print("Post Projection of particle1")
    #         print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
    #         print("Post Projection of particle2")
    #         print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
    #         print("")
    #         # print("Total Post Projection")
    #         # print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])+np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz]))
    #         # print("")
    #         print("Difference between pre and post")
    #         print("Particle 1")
    #         print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
    #         print("Particle 2")
    #         print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
    #         print("")

    #         gs2simdatalist[step2-1] = gausss2(startvec=postcollision,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_change)
    #         step2 -= 1
    #     # If the collision happens to be within tolerance
    #     else:
    #         print("Collision occurred within tolerance")

    #         prepart1x = s3kvecproj(gs2simdatalist[step2-1][0:3],gs2simdatalist[step2-1][6:9],"x")
    #         prepart1y = s3kvecproj(gs2simdatalist[step2-1][0:3],gs2simdatalist[step2-1][6:9],"y")
    #         prepart1z = s3kvecproj(gs2simdatalist[step2-1][0:3],gs2simdatalist[step2-1][6:9],"z")

    #         prepart2x = s3kvecproj(gs2simdatalist[step2-1][3:6],gs2simdatalist[step2-1][9:12],"x")
    #         prepart2y = s3kvecproj(gs2simdatalist[step2-1][3:6],gs2simdatalist[step2-1][9:12],"y")
    #         prepart2z = s3kvecproj(gs2simdatalist[step2-1][3:6],gs2simdatalist[step2-1][9:12],"z")

    #         print("Pre Projection of particle1")
    #         print([prepart1x,prepart1y,prepart1z])
    #         print("Pre Projection of particle2")
    #         print([prepart2x,prepart2y,prepart2z])
    #         print("")
    #         print("Total Pre Projection")
    #         print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
    #         print("")

    #         postcollision = h3collision(gs2simdatalist[step2-1],params)

    #         postpart1x = s3kvecproj(postcollision[0:3],postcollision[6:9],"x")
    #         postpart1y = s3kvecproj(postcollision[0:3],postcollision[6:9],"y")
    #         postpart1z = s3kvecproj(postcollision[0:3],postcollision[6:9],"z")

    #         postpart2x = s3kvecproj(postcollision[3:6],postcollision[9:12],"x")
    #         postpart2y = s3kvecproj(postcollision[3:6],postcollision[9:12],"y")
    #         postpart2z = s3kvecproj(postcollision[3:6],postcollision[9:12],"z")

    #         print("Post Projection of particle1")
    #         print([postpart1x,postpart1y,postpart1z])
    #         print("Post Projection of particle2")
    #         print([postpart2x,postpart2y,postpart2z])
    #         print("")
    #         print("Total Post Projection")
    #         print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z]))
    #         print("")
    #         print("Difference between pre and post")
    #         print(np.array([postpart1x,postpart1y,postpart1z])+np.array([postpart2x,postpart2y,postpart2z])-np.array([prepart1x,prepart1y,prepart1z])-np.array([prepart2x,prepart2y,prepart2z]))
    #         print("")

    #         gs2simdatalist[step2] = gausss2(startvec=postcollision,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)
    # else:
    #     # print("Not collided")
    #     gs2simdatalist[step2] = gausss2(startvec=gs2simdatalist[step2-1],params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)







    if abs((distck3 - (r1 +r2))) < tol or (distck3 - (r1 +r2)) < 0:
        print("Collided s3")
        dt_check = dt/2.
        # tol = 1e-6
        nump = 0

        if abs((distck3 - (r1 +r2))) > tol:
            print("Needed to find the collision position")
        
            dt_change = dt_check
            while abs((distck3 - (r1 +r2))) > tol and nump < 100:
                precollision = gausss3(startvec=gs3simdatalist[step3-2],params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt_change)
                distck3 = r4dist(rot2r4(precollision[0:3]),rot2r4(precollision[3:6]))
                if (distck3 - (r1 +r2)) < 0:
                    dt_change = dt_change - dt_check
                    # print("dt_check minus")
                    # print(dt_check)
                elif (distck3 - (r1 +r2)) > 0:
                    dt_check = dt_check/2
                    dt_change = dt_change + dt_check
                    # print("dt_check plus")
                    # print(dt_check)
                nump += 1

            prepart1bx, prepart1bxvec = s3kvecproj(precollision[0:3],precollision[6:9],"bx")
            prepart1by, prepart1byvec = s3kvecproj(precollision[0:3],precollision[6:9],"by")
            prepart1bz, prepart1bzvec = s3kvecproj(precollision[0:3],precollision[6:9],"bz")
            prepart1vx, prepart1vxvec = s3kvecproj(precollision[0:3],precollision[6:9],"vx")
            prepart1vy, prepart1vyvec = s3kvecproj(precollision[0:3],precollision[6:9],"vy")
            prepart1vz, prepart1vzvec = s3kvecproj(precollision[0:3],precollision[6:9],"vz")

            prepart2bx, prepart2bxvec = s3kvecproj(precollision[3:6],precollision[9:12],"bx")
            prepart2by, prepart2byvec = s3kvecproj(precollision[3:6],precollision[9:12],"by")
            prepart2bz, prepart2bzvec = s3kvecproj(precollision[3:6],precollision[9:12],"bz")
            prepart2vx, prepart2vxvec = s3kvecproj(precollision[3:6],precollision[9:12],"vx")
            prepart2vy, prepart2vyvec = s3kvecproj(precollision[3:6],precollision[9:12],"vy")
            prepart2vz, prepart2vzvec = s3kvecproj(precollision[3:6],precollision[9:12],"vz")

            print("Pre Projection of particle1")
            print([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])
            print("Pre Projection of particle2")
            print([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz])
            print("")
            # print("Total Pre Projection")
            # print(np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz])+np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
            # print("")

            postcollision = s3collision(precollision,params)

            postpart1bx, postpart1bxvec = s3kvecproj(postcollision[0:3],postcollision[6:9],"bx")
            postpart1by, postpart1byvec = s3kvecproj(postcollision[0:3],postcollision[6:9],"by")
            postpart1bz, postpart1bzvec = s3kvecproj(postcollision[0:3],postcollision[6:9],"bz")
            postpart1vx, postpart1vxvec = s3kvecproj(postcollision[0:3],postcollision[6:9],"vx")
            postpart1vy, postpart1vyvec = s3kvecproj(postcollision[0:3],postcollision[6:9],"vy")
            postpart1vz, postpart1vzvec = s3kvecproj(postcollision[0:3],postcollision[6:9],"vz")

            postpart2bx, postpart2bxvec = s3kvecproj(postcollision[3:6],postcollision[9:12],"bx")
            postpart2by, postpart2byvec = s3kvecproj(postcollision[3:6],postcollision[9:12],"by")
            postpart2bz, postpart2bzvec = s3kvecproj(postcollision[3:6],postcollision[9:12],"bz")
            postpart2vx, postpart2vxvec = s3kvecproj(postcollision[3:6],postcollision[9:12],"vx")
            postpart2vy, postpart2vyvec = s3kvecproj(postcollision[3:6],postcollision[9:12],"vy")
            postpart2vz, postpart2vzvec = s3kvecproj(postcollision[3:6],postcollision[9:12],"vz")

            print("Post Projection of particle1")
            print([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])
            print("Post Projection of particle2")
            print([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])
            print("")
            # print("Total Post Projection")
            # print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])+np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz]))
            # print("")
            print("Difference between pre and post")
            print("Particle 1")
            print(np.array([postpart1bx,postpart1by,postpart1bz,postpart1vx,postpart1vy,postpart1vz])-np.array([prepart1bx,prepart1by,prepart1bz,prepart1vx,prepart1vy,prepart1vz]))
            print("Particle 2")
            print(np.array([postpart2bx,postpart2by,postpart2bz,postpart2vx,postpart2vy,postpart2vz])-np.array([prepart2bx,prepart2by,prepart2bz,prepart2vx,prepart2vy,prepart2vz]))
            print("")

            gs3simdatalist[step3-1] = gausss3(startvec=postcollision,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_h3sim2ballcol,dt=dt-dt_change)
            step3 -= 1
        # If the collision happens to be within tolerance
        else:
            print("Collision occurred within tolerance")

            prepart1x = s3kvecproj(gs3simdatalist[step3-1][0:3],gs3simdatalist[step3-1][6:9],"vx")
            prepart1y = s3kvecproj(gs3simdatalist[step3-1][0:3],gs3simdatalist[step3-1][6:9],"vy")
            prepart1z = s3kvecproj(gs3simdatalist[step3-1][0:3],gs3simdatalist[step3-1][6:9],"vz")

            prepart2x = s3kvecproj(gs3simdatalist[step3-1][3:6],gs3simdatalist[step3-1][9:12],"vx")
            prepart2y = s3kvecproj(gs3simdatalist[step3-1][3:6],gs3simdatalist[step3-1][9:12],"vy")
            prepart2z = s3kvecproj(gs3simdatalist[step3-1][3:6],gs3simdatalist[step3-1][9:12],"vz")

            print("Pre Projection of particle1")
            print([prepart1x,prepart1y,prepart1z])
            print("Pre Projection of particle2")
            print([prepart2x,prepart2y,prepart2z])
            print("")
            print("Total Pre Projection")
            print(np.array([prepart1x,prepart1y,prepart1z])+np.array([prepart2x,prepart2y,prepart2z]))
            print("")

            postcollision = h3collision(gs3simdatalist[step3-1],params)

            postpart1x = s3kvecproj(postcollision[0:3],postcollision[6:9],"vx")
            postpart1y = s3kvecproj(postcollision[0:3],postcollision[6:9],"vy")
            postpart1z = s3kvecproj(postcollision[0:3],postcollision[6:9],"vz")

            postpart2x = s3kvecproj(postcollision[3:6],postcollision[9:12],"vx")
            postpart2y = s3kvecproj(postcollision[3:6],postcollision[9:12],"vy")
            postpart2z = s3kvecproj(postcollision[3:6],postcollision[9:12],"vz")

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

            gs3simdatalist[step3] = gausss2(startvec=postcollision,params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)
    else:
        # print("Not collided")
        gs3simdatalist[step3] = gausss3(startvec=gs3simdatalist[step3-1],params=params,dynfunc=dynfunc_s3sim2ballcol,dynjac=dynjac_s3sim2ballcol,dt=dt)


    # Sim in H3
    # # startvecgs1sim = gs1simdatalist[step1]
    # # startvecgs2sim = gs2simdatalist[step2]
    # startvecgs3sim = gs3simdatalist[step3]

    # # distck1 = h3dist(rot2hyp(startvecgs1sim[0:3]),rot2hyp(startvecgs1sim[3:6]))
    # # distck2 = h3dist(rot2hyp(startvecgs2sim[0:3]),rot2hyp(startvecgs2sim[3:6]))
    # distck3 = h3dist(rot2hyp(startvecgs3sim[0:3]),rot2hyp(startvecgs3sim[3:6]))

    # Sim in S3
    # startvecgs1sim = gs1simdatalist[step1]
    # startvecgs2sim = gs2simdatalist[step2]
    startvecgs3sim = gs3simdatalist[step3]

    # distck1 = r4dist(rot2r4(gs1simdatalist[step1][0:3]),rot2r4(gs1simdatalist[step1][3:6]))
    # distck2 = r4dist(rot2r4(gs2simdatalist[step2][0:3]),rot2r4(gs2simdatalist[step2][3:6]))
    distck3 = r4dist(rot2r4(gs3simdatalist[step3][0:3]),rot2r4(gs3simdatalist[step3][3:6]))

    # print("dist3")
    # print(distck3)
    # print("dist3 - r1 -r2")
    # print((distck3 - (r1 +r2)))


    # print(step3)
    # if step1%int(1/dt)==0:
    #         print(step1)
    # if step2%int(1/dt)==0:
    #         print(step2)
    if step3%int(1/dt)==0:
            print(step3)
    # step1 += 1
    # step2 += 1
    step3 += 1


# Sim in H3
# np.save("h3_1d_col_gausss1_sim_tmax10_dt1",gs1simdatalist)
# np.save("h3_1d_col_gausss2_sim_tmax10_dt1",gs2simdatalist)
np.save("h3_1d_col_gausss3_sim_tmax10_dt1",gs3simdatalist)

# data1 = np.load("h3_1d_col_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("h3_1d_col_gausss2_sim_tmax10_dt1.npy")
data3 = np.load("h3_1d_col_gausss3_sim_tmax10_dt1.npy")

# np.save("h3_2d_col_gausss1_sim_tmax10_dt1",gs1simdatalist)
# np.save("h3_2d_col_gausss2_sim_tmax10_dt1",gs2simdatalist)
# np.save("h3_2d_col_gausss3_sim_tmax10_dt1",gs3simdatalist)

# data1 = np.load("h3_2d_col_gausss1_sim_tmax10_dt1.npy")
# data2 = np.load("h3_2d_col_gausss2_sim_tmax10_dt1.npy")
# data3 = np.load("h3_2d_col_gausss3_sim_tmax10_dt1.npy")

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
energydata3 = np.zeros(np.shape(data3)[0])

counter = 0
for a in range(np.shape(data3)[0]):
#     distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
#     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
    distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
    energydata3[counter] = .5*m1*(data3[a][6]**2.+np.sin(data3[a][0])**2.*data3[a][7]**2.+np.sin(data3[a][0])**2.*np.sin(data3[a][1])**2.*data3[a][8]**2.) +.5*m2*(data3[a][9]**2.+np.sin(data3[a][3])**2.*data3[a][10]**2.+np.sin(data3[a][3])**2.*np.sin(data3[a][4])**2.*data3[a][11]**2.)
    counter += 1

# counter = 0
# for a in range(np.shape(data3)[0]):
#     # distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
#     # distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
#     distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
#     energydata3[counter] = .5*m1*(data3[a][6]**2.+np.sinh(data3[a][0])**2.*data3[a][7]**2.+np.sinh(data3[a][0])**2.*np.sin(data3[a][1])**2.*data3[a][8]**2.) +.5*m2*(data3[a][9]**2.+np.sinh(data3[a][3])**2.*data3[a][10]**2.+np.sinh(data3[a][3])**2.*np.sin(data3[a][4])**2.*data3[a][11]**2.)
#     counter += 1

# ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")
# ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# ax.plot(t_arr,distdata1,'r',label = "Gauss s1")
# ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
ax.plot(t_arr,np.full(len(t_arr),1.))
# ax.plot(t_arr[:counter-1],energydata3[:counter-1])
ax.legend()
ax.set_title('Simulation Data')
# ax.set_xlim(2.95,3.05)
# ax.set_ylim(.95,1.05)
ax.set_xlabel('t')
ax.set_ylabel('l')
plt.show()







