import numpy as np
# from rod_derivative_bank import *
from src.integrator_files.integrator_bank import h3rads2roddae, h3rads3roddae, s3rads2roddae, s3rads3roddae
from src.test_system_simulations.test_system_bank import dynfunc_h3simbar, dynjac_h3simbar, dynfunc_s3simbar, dynjac_s3simbar



# Simulation class setup

class RigidBarSimulation:

    def __init__(self, ambient_geo, system_params, dt, tmax, solver_id):
        # Ambient Space
        self.ambient_geo = ambient_geo

        # Time Data
        self.dt = dt   
        self.tmax = tmax
        self.t_arr = np.arange(0.,self.tmax+self.dt,self.dt)

        # Test System Data
        self.simdatalist = np.zeros((self.t_arr.shape[0],12+1))

        # System Parameters [ m1, m2, x ]
        self.system_params  = system_params

        # Integrator to use
        self.solver_id = solver_id
        self.tol = 1e-10

        # Internal Flags
        self.have_ics = False
        self.have_run = False

    def set_initial_conditions(self, system_ics):
        self.system_ics = system_ics
        self.have_ics = True

    def clear_data(self):
        self.simdatalist = np.zeros((self.t_arr.shape[0],12+1))
        self.have_run = False

    def run(self):
        if self.have_ics:
            # Hyperbolic Space (Sim)
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                self.simdatalist[0] = self.system_ics.copy()
                        
                # Radau 2-Step Method (Based off the work of Schweizer and Li)
                if self.solver_id == "rs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = h3rads2roddae(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 3-Step Method (Based off the work of Schweizer and Li)
                if self.solver_id == "rs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = h3rads3roddae(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dt=self.dt,
                            tol=self.tol)

            # Spherical Space (Sim)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                self.simdatalist[0] = self.system_ics.copy()

                        
                # Radau 2-Step Method (Based off the work of Schweizer and Li)
                if self.solver_id == "rs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = s3rads2roddae(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 3-Step Method (Based off the work of Schweizer and Li)
                if self.solver_id == "rs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = s3rads3roddae(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dt=self.dt,
                            tol=self.tol)
                        
            print("Simulation run completed!")
            self.have_run = True
        else:
            print("Error: Must provide initial conditions via set_initial_conditions() before running simulation")

    def output_data(self):
        if self.have_run:
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                np.save("h3_r_rig_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                np.save("s3_r_rig_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
        else:
            print("Error: Must use run() to generate data")



# import matplotlib.pyplot as plt
# from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3

# # Solver Setup

# # Time array based on time step
# dt = .1    # Number of steps
# t_max = 10.      # Total simulation time

# # Initial Data
# v = 1.      # Initial Velocity
# # x = 1.      # Spring Rest Length H3
# x = 1.      # Rod Length H3
# m1 = 1.      # Mass of point masses
# m2 = 1.      # Mass of point masses
# params = [m1,m2,x]


# # Sim bar in H3
# geometry = "h3"
# startvec = np.array([
#     [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()
# startvec = np.append(startvec,0.5876005968219006)
# # startvec = np.append(startvec,0.) # lambda=0 old guess

# # Sim bar in S3
# # startvec = np.array([
# #     [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
# #     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()
# # startvec = np.append(startvec,-0.4207354924039482)
# # startvec = np.append(startvec,0.) # lambda=0 old guess

# # Solver 
# solver_id = "rs2"

# # Initialize Simulation Object
# sim_test = RigidBarSimulation(geometry, params, dt, t_max, solver_id, True)
# sim_test.set_initial_conditions(startvec)
# sim_test.run()
# sim_test.output_data()

# data1 = np.load("h3_r_rig_rs2_sim_tmax10.0_dt0.1.npy")


# fig,ax=plt.subplots(1,1)

# distdata1 = np.zeros(np.shape(data1)[0])

# # counter = 0
# # for a in range(np.shape(data2)[0]):
# #     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))-1.
# #     distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))-1.
# #     counter += 1

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
#     # distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
#     counter += 1

# # ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# # ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# # ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")
# # ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# # ax.plot(t_arr,distdata1,'r',label = "Gauss s1")

# ax.plot(sim_test.t_arr,distdata1,marker = ".",color='k',linestyle = "None",label = "Radau s2")
# # ax.plot(t_arr,distdata3,marker = ".",color='b',linestyle = "None",label = "Radau s3")
# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()













        
    