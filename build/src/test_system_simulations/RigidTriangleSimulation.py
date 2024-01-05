import numpy as np
from build.src.integrator_files.integrator_bank import h3rads2tridae, h3rads3tridae, s3rads2tridae, s3rads3tridae
from build.src.test_system_simulations.test_system_bank import dynfunc_h3simtriangle, dynjac_h3simtriangle, dynfunc_s3simtriangle, dynjac_s3simtriangle


# Simulation class setup

class RigidTriangleSimulation:

    def __init__(self, ambient_geo, system_params, dt, tmax, solver_id):
        # Ambient Space
        self.ambient_geo = ambient_geo

        # Time Data
        self.dt = dt   
        self.tmax = tmax
        self.t_arr = np.arange(0.,self.tmax+self.dt,self.dt)

        # Test System Data
        self.simdatalist = np.zeros((self.t_arr.shape[0],18+3))

        # System Parameters [ m1, m2, m3, x1, x2, x3 ]
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
        self.simdatalist = np.zeros((self.t_arr.shape[0],18+3))
        self.have_run = False

    def run(self):
        if self.have_ics:
            # Hyperbolic Space (Sim)
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                self.simdatalist[0] = self.system_ics.copy()
                        
                # Radau 2-Step Method (Based off the work of Schweizer and Li)
                if self.solver_id == "rs2":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = h3rads2tridae(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 3-Step Method (Based off the work of Schweizer and Li)
                if self.solver_id == "rs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = h3rads3tridae(
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
                        self.simdatalist[step+1] = s3rads2tridae(
                            startvec=self.simdatalist[step],
                            params=self.system_params,
                            dt=self.dt,
                            tol=self.tol)
                        
                # Radau 3-Step Method (Based off the work of Schweizer and Li)
                if self.solver_id == "rs3":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.simdatalist[step+1] = s3rads3tridae(
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
                np.save("h3_t_rig_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                np.save("s3_t_rig_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
        else:
            print("Error: Must use run() to generate data")




# import matplotlib.pyplot as plt
# from function_bank import rot2hyp, hyp2rot, hyp2poin3d, h3dist, killingvech3, rot2r4, r42rot, s2rstproj, r4dist, killingvecs3

# dt = .1         # Number of steps
# t_max = 10.      # Total simulation time

# # Initial Data
# v = 1.      # Initial Velocity
# # x = 1.      # Spring Rest Length H3
# x1 = 1.      # Rod Length H3
# x2 = 1.      # Rod Length H3
# x3 = 1.      # Rod Length H3
# m1 = 1.      # Mass of point masses
# m2 = 1.      # Mass of point masses
# m3 = 1.      # Mass of point masses
# params = [m1,m2,m3,x1,x2,x3]


# # Sim bar in H3
# geometry = "h3"
# startvec = np.array([
#     [np.arccosh(np.cosh(x1)/np.cosh(x1/2.)),np.pi/2.,0.],[.5,np.pi/2.,np.pi/2.],[.5,np.pi/2.,3.*np.pi/2.],
#     killingvech3([np.arccosh(np.cosh(x1)/np.cosh(x1/2.)),np.pi/2.,0.],v,"x"), killingvech3([.5,np.pi/2.,np.pi/2.],v,"x"), killingvech3([.5,np.pi/2.,3.*np.pi/2.],v,"x")]).flatten()
# startvec = np.append(startvec,0.) #lam1
# startvec = np.append(startvec,0.) #lam2
# startvec = np.append(startvec,0.5876005968219006)
# # # startvec = np.append(startvec,0.) #lam3  initial zero

# # Sim bar in S3
# # startvec = np.array([
# #     [np.pi/2., np.pi/2., np.arccos(cos(x1)/cos(x1/2.))],[(np.pi - 1.)/2.,np.pi/2.,0.],[(np.pi + 1.)/2.,np.pi/2.,0.],
# #     killingvecs3([np.pi/2., np.pi/2., np.arccos(cos(x1)/cos(x1/2.))],-v,"vz"), killingvecs3([(np.pi - 1.)/2.,np.pi/2.,0.],-v,"vz"), killingvecs3([(np.pi + 1.)/2.,np.pi/2.,0.],-v,"vz")]).flatten()
# # startvec = np.append(startvec,0.) #lam1
# # startvec = np.append(startvec,0.) #lam2
# # startvec = np.append(startvec,-0.4207354924039482)
# # startvec = np.append(startvec,0.) #lam3 initial zero

# # Solver 
# solver_id = "rs2"

# # Initialize Simulation Object
# sim_test = RigidTriangleSimulation(geometry, params, dt, t_max, solver_id)
# sim_test.set_initial_conditions(startvec)
# sim_test.run()
# sim_test.output_data()

# data1 = np.load("h3_t_rig_rs2_sim_tmax10.0_dt0.1.npy")

# fig,ax=plt.subplots(1,1)

# distdatasp12r2 = np.zeros(np.shape(data1)[0])
# distdatasp13r2 = np.zeros(np.shape(data1)[0])
# distdatasp23r2 = np.zeros(np.shape(data1)[0])
# # distdatasp12r3 = np.zeros(np.shape(data1)[0])
# # distdatasp13r3 = np.zeros(np.shape(data1)[0])
# # distdatasp23r3 = np.zeros(np.shape(data1)[0])

# counter = 0
# for a in range(np.shape(data1)[0]):
#     distdatasp12r2[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][3:6]))-1.
#     distdatasp13r2[counter] = h3dist(rot2hyp(data1[a][0:3]),rot2hyp(data1[a][6:9]))-1.
#     distdatasp23r2[counter] = h3dist(rot2hyp(data1[a][3:6]),rot2hyp(data1[a][6:9]))-1.
#     # distdatasp13r3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][6:9]))-1.
#     # distdatasp12r3[counter] = h3dist(rot2hyp(data3[a][0:3]),rot2hyp(data3[a][3:6]))-1.
#     # distdatasp23r3[counter] = h3dist(rot2hyp(data3[a][3:6]),rot2hyp(data3[a][6:9]))-1.
#     counter += 1


# # ax.plot(t_arr,2.*(np.pi/2. - gs1exactdatalist[:,0]),'r',label = "Gauss s1")
# # ax.plot(t_arr,2.*(np.pi/2. - gs2exactdatalist[:,0]),'k',label = "Gauss s2")
# # ax.plot(t_arr,2.*(np.pi/2. - gs3exactdatalist[:,0]),'b',label = "Gauss s3")
# # ax.plot(t_arr,2.*(data1[:,0]),'b',label = "Gauss h3")
# # ax.plot(t_arr,distdata1,'r',label = "Gauss s1")

# ax.plot(sim_test.t_arr,distdatasp12r2,marker = ".",color='k',linestyle = "None",label = "Radau s2 d12")
# ax.plot(sim_test.t_arr,distdatasp13r2,marker = ".",color='k',linestyle = "None",label = "Radau s2 d13")
# ax.plot(sim_test.t_arr,distdatasp23r2,marker = ".",color='k',linestyle = "None",label = "Radau s2 d23")
# # ax.plot(t_arr,distdatasp12r3,marker = ".",color='b',linestyle = "None",label = "Radau s3 d12")
# # ax.plot(t_arr,distdatasp13r3,marker = ".",color='b',linestyle = "None",label = "Radau s3 d13")
# # ax.plot(t_arr,distdatasp23r3,marker = ".",color='b',linestyle = "None",label = "Radau s3 d23")

# ax.legend()
# ax.set_title('Simulation Data')
# ax.set_xlabel('t')
# ax.set_ylabel('l')
# plt.show()
    