import numpy as np
from src.integrator_files.integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from src.test_system_simulations.test_system_bank import dynfunc_h3simbar, dynjac_h3simbar, dynfunc_s3simbar, dynjac_s3simbar


# Simulation class setup

class SpringBarSimulation:

    def __init__(self, ambient_geo, system_params, dt, tmax, solver_id):
        # Ambient Space
        self.ambient_geo = ambient_geo

        # Time Data
        self.dt = dt   
        self.tmax = tmax
        self.t_arr = np.arange(0.,self.tmax+self.dt,self.dt)

        # Test System Data
        self.simdatalist = np.zeros((self.t_arr.shape[0],12))

        # System Parameters [ v, ks, x, m ]
        self.system_params  = system_params

        # Integrator to use
        self.solver_id = solver_id
        self.tol = 1e-15

        # Internal Flags
        self.have_ics = False
        self.have_run = False

    def set_initial_conditions(self, system_ics):
        self.system_ics = system_ics
        self.have_ics = True

    def clear_data(self):
        self.simdatalist = np.zeros((self.t_arr.shape[0],12))
        self.have_run = False

    def run(self):
        if self.have_ics:
            # Hyperbolic Space (Sim)
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                self.simdatalist[0] = self.system_ics.copy()

                # Gauss 1-Step Method
                if self.solver_id == "gs1":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss1(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simbar,
                            dynjac=dynjac_h3simbar,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 2-Step Method
                if self.solver_id == "gs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss2(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simbar,
                            dynjac=dynjac_h3simbar,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 3-Step Method
                if self.solver_id == "gs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss3(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simbar,
                            dynjac=dynjac_h3simbar,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 2-Step Method
                if self.solver_id == "rs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = rads2(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simbar,
                            dynjac=dynjac_h3simbar,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 3-Step Method
                if self.solver_id == "rs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = rads3(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simbar,
                            dynjac=dynjac_h3simbar,
                            dt=self.dt,
                            tol=self.tol)

            # Spherical Space (Sim)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                self.simdatalist[0] = self.system_ics.copy()

                # Gauss 1-Step Method
                if self.solver_id == "gs1":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss1(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_s3simbar,
                            dynjac=dynjac_s3simbar,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 2-Step Method
                if self.solver_id == "gs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss2(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_s3simbar,
                            dynjac=dynjac_s3simbar,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 3-Step Method
                if self.solver_id == "gs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss3(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_s3simbar,
                            dynjac=dynjac_s3simbar,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 2-Step Method
                    if self.solver_id == "rs2":
                        for step in range(self.t_arr.shape[0] - 1):
                            self.simdatalist[step+1] = rads2(
                                startvec=self.simdatalist[step],
                                params=self.system_params,
                                dynfunc=dynfunc_s3simbar,
                                dynjac=dynjac_s3simbar,
                                dt=self.dt,
                                tol=self.tol)
                            
                    # Radau 3-Step Method
                    if self.solver_id == "rs3":
                        for step in range(self.t_arr.shape[0] - 1):
                            self.simdatalist[step+1] = rads3(
                                startvec=self.simdatalist[step],
                                params=self.system_params,
                                dynfunc=dynfunc_s3simbar,
                                dynjac=dynjac_s3simbar,
                                dt=self.dt,
                                tol=self.tol)

            print("Simulation run completed!")
            self.have_run = True
        else:
            print("Error: Must provide initial conditions via set_initial_conditions() before running simulation")

    def output_data(self):
        if self.have_run:
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                np.save("h3_r_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                np.save("s3_r_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
        else:
            print("Error: Must use run() to generate data")

        

# import matplotlib.pyplot as plt
# from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3

# # Initialize Test System
                
# dt = .1         # Number of steps
# t_max = 10.      # Total simulation time

# # Initial Data
# v = 1.      # Initial Velocity
# ks = 1.     # Spring Stiffness
# x = 1.      # Spring Rest Length H3 and S3
# # x = np.pi - 1.      # Spring Rest Length S3 exact
# m = 1.      # Mass of point masses
# params = [v,ks,x,m]


# # Exact bar in H3
# # geometry, system_id = "h3", "exact"
# # startvec = np.array([.5,0.])
# # Sim bar in H3
# geometry = "h3"
# startvec = np.array([
#     [.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()

# # Exact bar in S3
# # geometry, system_id = "s3", "exact"
# # startvec = np.array([(np.pi - 1.)/2.,0.])
# # Sim bar in S3
# # geometry, system_id = "s3", "sim"
# # startvec = np.array([
# #     [(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
# #     killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()

# # Solver 
# solver_id = "gs1"

# # Initialize Simulation Object
# sim_test = SpringBarSimulation(geometry, params, dt, t_max, solver_id)
# sim_test.set_initial_conditions(startvec)
# sim_test.run()
# sim_test.output_data()

# data1 = np.load("h3_r_gs1_sim_tmax10.0_dt0.1.npy")

# fig,ax=plt.subplots(1,1)

# distdata1 = np.zeros(np.shape(data1)[0])
# # distdata2 = np.zeros(np.shape(data2)[0])
# # distdata3 = np.zeros(np.shape(data3)[0])

# # counter = 0
# # for a in range(np.shape(data1)[0]):
# #     distdata1[counter] = r4dist(rot2r4(data1[a][0:3]),rot2r4(data1[a][3:6]))
# #     distdata2[counter] = r4dist(rot2r4(data2[a][0:3]),rot2r4(data2[a][3:6]))
# #     distdata3[counter] = r4dist(rot2r4(data3[a][0:3]),rot2r4(data3[a][3:6]))
# #     counter += 1

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdata1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
#     # distdata2[counter] = h3dist(rot2hyp(data2[a][0:3]),rot2hyp(data2[a][3:6]))
#     # distdata3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))
#     counter += 1

# # ax.plot(t_arr,2.*(np.pi/2. - data1[:,0]),'r',label = "Gauss s1")
# # ax.plot(t_arr,2.*(np.pi/2. - data2[:,0]),'k',label = "Gauss s2")
# # ax.plot(t_arr,2.*(np.pi/2. - data3[:,0]),'b',label = "Gauss s3")
# # ax.plot(sim_test.t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# # ax.plot(t_arr,2.*(data2[:,0]),'b',label = "Gauss h3")
# # ax.plot(t_arr,2.*(data3[:,0]),'b',label = "Gauss h3")
# ax.plot(sim_test.t_arr,distdata1,'r',label = "Gauss s1")
# # ax.plot(t_arr,distdata2,'k',label = "Gauss s2")
# # ax.plot(t_arr,distdata3,'b',label = "Gauss s3")
# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()
    