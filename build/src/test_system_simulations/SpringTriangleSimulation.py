import numpy as np
from build.src.integrator_files.integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from build.src.test_system_simulations.test_system_bank import dynfunc_h3simtriangle, dynjac_h3simtriangle, dynfunc_s3simtriangle, dynjac_s3simtriangle


# Simulation class setup

class SpringTriangleSimulation:

    def __init__(self, ambient_geo, system_params, dt, tmax, solver_id):
        # Ambient Space
        self.ambient_geo = ambient_geo

        # Time Data
        self.dt = dt   
        self.tmax = tmax
        self.t_arr = np.arange(0.,self.tmax+self.dt,self.dt)

        # Test System Data
        self.simdatalist = np.zeros((self.t_arr.shape[0],18))

        # System Parameters [ v, ks, x, m ]
        self.system_params  = system_params

        # Integrator to use
        self.solver_id = solver_id
        self.tol = 1e-15

        self.have_ics = False
        self.have_run = False

    def set_initial_conditions(self, system_ics):
        self.system_ics = system_ics
        self.have_ics = True

    def clear_data(self):
        self.simdatalist = np.zeros((self.t_arr.shape[0],18))
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
                            dynfunc=dynfunc_h3simtriangle,
                            dynjac=dynjac_h3simtriangle,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 2-Step Method
                if self.solver_id == "gs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss2(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simtriangle,
                            dynjac=dynjac_h3simtriangle,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 3-Step Method
                if self.solver_id == "gs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss3(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simtriangle,
                            dynjac=dynjac_h3simtriangle,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 2-Step Method
                if self.solver_id == "rs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = rads2(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simtriangle,
                            dynjac=dynjac_h3simtriangle,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 3-Step Method
                if self.solver_id == "rs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = rads3(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_h3simtriangle,
                            dynjac=dynjac_h3simtriangle,
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
                            dynfunc=dynfunc_s3simtriangle,
                            dynjac=dynjac_s3simtriangle,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 2-Step Method
                if self.solver_id == "gs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss2(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_s3simtriangle,
                            dynjac=dynjac_s3simtriangle,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Gauss 3-Step Method
                if self.solver_id == "gs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = gausss3(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dynfunc=dynfunc_s3simtriangle,
                            dynjac=dynjac_s3simtriangle,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 2-Step Method
                    if self.solver_id == "rs2":
                        for step in range(self.t_arr.shape[0] - 1):
                            self.simdatalist[step+1] = rads2(
                                startvec=self.simdatalist[step],
                                params=self.system_params,
                                dynfunc=dynfunc_s3simtriangle,
                                dynjac=dynjac_s3simtriangle,
                                dt=self.dt,
                                tol=self.tol)
                            
                    # Radau 3-Step Method
                    if self.solver_id == "rs3":
                        for step in range(self.t_arr.shape[0] - 1):
                            self.simdatalist[step+1] = rads3(
                                startvec=self.simdatalist[step],
                                params=self.system_params,
                                dynfunc=dynfunc_s3simtriangle,
                                dynjac=dynjac_s3simtriangle,
                                dt=self.dt,
                                tol=self.tol)

            print("Simulation run completed!")
            self.have_run = True
        else:
            print("Error: Must provide initial conditions via set_initial_conditions() before running simulation")

    def output_data(self):
        if self.have_run:
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                np.save("h3_t_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                np.save("s3_t_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
        else:
            print("Error: Must use run() to generate data")



# import matplotlib.pyplot as plt
# from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3

# dt = .1         # Number of steps
# t_max = 10.      # Total simulation time

# # Initial Data
# v = 1.      # Initial Velocity
# ks = 1.     # Spring Stiffness
# # x = 1.      # Spring Rest Length H3
# # x = np.pi - 1.      # Spring Rest Length S3 exact
# x = 1.      # Spring Rest Length S3
# m = 1.      # Mass of point masses
# params = [v,ks,x,m]

# # Exact triangle in H3
# # startvec = np.array([x/2.,np.arccosh(np.cosh(x)/np.cosh(x/2.)),0.,0.])
# # Sim triangle in H3
# geometry = "h3"
# startvec = np.array([
#     [np.arccosh(np.cosh(x)/np.cosh(x/2.)),np.pi/2.,0.],[.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([np.arccosh(np.cosh(x)/np.cosh(x/2.)),np.pi/2.,0.],v,"x"), killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()
# # Exact triangle in S3
# # startvec = np.array([x/2.,np.arccos(np.cos(x)/np.cos(x/2.)),0.,0.])
# # Sim bar in S3
# # startvec = np.array([
# #     [np.pi/2.,np.pi/2.,np.arccos(np.cos(x)/np.cos(x/2.))],[(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0],
# #     killingvecs3([np.pi/2.,np.pi/2.,np.arccos(np.cos(x)/np.cos(x/2.))],-v,"vz"), killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0],-v,"vz")]).flatten()

# # Solver 
# solver_id = "gs1"

# # Initialize Simulation Object
# sim_test = SpringTriangleSimulation(geometry, params, dt, t_max, solver_id)
# sim_test.set_initial_conditions(startvec)
# sim_test.run()
# sim_test.output_data()

# data1 = np.load("h3_t_gs1_sim_tmax10.0_dt0.1.npy")

# fig,ax=plt.subplots(1,1)

# distdatasp12g1 = np.zeros(np.shape(data1)[0])
# distdatasp13g1 = np.zeros(np.shape(data1)[0])
# distdatasp23g1 = np.zeros(np.shape(data1)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdatasp12g1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))
#     distdatasp13g1[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][6:9]))
#     distdatasp23g1[counter] = h3dist(rot2hyp(data1[a][3:6]),rot2hyp(data1[a][6:9]))
#     counter += 1


# ax.plot(sim_test.t_arr,distdatasp12g1,'r',label = "Gauss s1")

# ax.plot(sim_test.t_arr,distdatasp13g1,'k',label = "Gauss s1")

# ax.plot(sim_test.t_arr,distdatasp23g1,'b',label = "Gauss s1")

# # ax.plot(t_arr,2*data3[:,0],'b',label = "Gauss s3")
# # ax.plot(t_arr,np.arccosh(np.cosh(data3[:,0])*np.cosh(data3[:,1])),'r',label = "Gauss s3")

# # ax.plot(t_arr,2*data3[:,0],'b',label = "Gauss s3")
# # ax.plot(t_arr,np.arccos(np.cos(data3[:,0])*np.cos(data3[:,1])),'r',label = "Gauss s3")
# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_ylim(0,2.5)
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()
    