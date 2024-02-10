import numpy as np
from src.integrator_files.integrator_bank import gausss1, gausss2, gausss3, rads2, rads3
from src.test_system_simulations.test_system_bank import dynfunc_h3sim2ballcol, dynjac_h3sim2ballcol
from src.collision_function_bank import h3collision, s3collision


# Simulation class setup

class VertexCollisionSimulation:
    """
    A class used to perform vertex collision simulation with two free vertices

    ...

    Attributes
    ----------
    ambient_geo : str
        the geometry of ambient space in simulation environment
        - currently supporting:
            * h3   = 3D hyperbolic space
            * s3   = 3D spherical space

    system_params : array
        the array of parameters describing system consisting of:
            * m1 = mass of vertex 1
            * m2 = mass of vertex 2
            * r1 = radius of pseudo rigid shell of vertex 1
            * r2 = radius of pseudo rigid shell of vertex 2 

    dt : float
        the simuation time step size
    tmax : float
        the total simulation time

    solver_id : str
        the solver to be used to evaluate dynamics
        - currently supporting:
            * gs1 = 1-step Gauss collocation 
            * gs2 = 2-step Gauss collocation 
            * gs3 = 3-step Gauss collocation 
            * rs2 = 2-step Radau collocation (RadauIIA)
            * rs3 = 3-step Radau collocation (RadauIIA)

    Methods
    -------
    set_initial_conditions(system_ics)
        Inputs user given initial conditions of system for simulation

    clear_data()
        Clears any simulation data stored in simulation object

    run()
        Runs simulation once given all necessary information

    output_data()
        Outputs simulation data to file with name:
            {self.ambient_geo}_r_{self.solver_id}_sim_tmax{self.tmax}_dt{self.dt}.npy
    """

    def __init__(self, ambient_geo, system_params, dt, tmax, solver_id):
        """
        Parameters
        ----------
        ambient_geo : str
            the geometry of ambient space in simulation environment
            - currently supporting:
                * h3 = 3D hyperbolic space
                * s3 = 3D spherical space

        system_params : array
            the array of parameters describing system consisting of:
                * m1 = mass of vertex 1
                * m2 = mass of vertex 2
                * r1 = radius of pseudo rigid shell of vertex 1
                * r2 = radius of pseudo rigid shell of vertex 2   

        dt : float
            the simuation time step size
        tmax : float
            the total simulation time

        solver_id : str
            the solver to be used to evaluate dynamics
            - currently supporting:
                * gs1 = 1-step Gauss collocation 
                * gs2 = 2-step Gauss collocation 
                * gs3 = 3-step Gauss collocation 
                * rs2 = 2-step Radau collocation (RadauIIA)
                * rs3 = 3-step Radau collocation (RadauIIA)

        """

        # Ambient Space
        self.ambient_geo = ambient_geo

        # System Parameters [ m1, m2, r1, r2 ]
        self.system_params  = system_params

        # Time Data
        self.dt = dt   
        self.tmax = tmax
        self.t_arr = np.arange(0.,self.tmax+self.dt,self.dt)

        # Test System Data
        self.simdatalist = np.zeros((self.t_arr.shape[0],12))

        # Integrator to use
        self.solver_id = solver_id
        self.tol = 1e-15
        self.collision_tol = 1e-8
        self.max_collision_check = 100

        # Internal Flags
        self._have_ics = False
        self._have_run = False

    def check_for_collision(self,step):

        self.m1 = self.system_params[0]
        self.m2 = self.system_params[1]
        self.r1 = self.system_params[2]
        self.r2 = self.system_params[3]

        # H3 Collision Methods
        def h3dist(v1,v2):
            x1,y1,z1,w1 = v1
            x2,y2,z2,w2 = v2
            return np.arccosh(-x1*x2 - y1*y2 - z1*z2 + w1*w2)
        
        def rot2hyp(vec):
            a,b,g = vec
            return np.array([
                np.sinh(a)*np.sin(b)*np.cos(g), 
                np.sinh(a)*np.sin(b)*np.sin(g), 
                np.sinh(a)*np.cos(b), 
                np.cosh(a)])

        def hyp2rot(vec): 
            x,y,z,w = vec
            return np.array([
                np.arccosh(w), 
                np.arccos(z/np.sinh(np.arccosh(w))), 
                np.arctan2(y, x)])
        
        # S3 Collision Methods
        def r4dist(v1, v2):
            x1,y1,z1,w1 = v1
            x2,y2,z2,w2 = v2
            return np.arccos(x1*x2 + y1*y2 + z1*z2 + w1*w2)
        
        def rot2r4(vec): 
            a,b,g = vec
            return np.array([
                np.sin(a)*np.sin(b)*np.cos(g),
                np.sin(a)*np.sin(b)*np.sin(g), 
                np.sin(a)*np.cos(b),
                np.cos(a)])

        def r42rot(vec): 
            x,y,z,w = vec
            return np.array([
                np.arccos(w), 
                np.arccos(z/np.sin(np.arccos(w))),
                np.arctan2(y, x)])
        

        # H3 collision check
        if self.ambient_geo == "h3" or self.ambient_geo == "H3":

            self.distance_check = h3dist(rot2hyp(self.simdatalist[step][0:3]),rot2hyp(self.simdatalist[step][3:6]))

            if abs((self.distance_check - (self.r1 + self.r2))) < self.collision_tol or self.distance_check - (self.r1 + self.r2) < 0.:
                print("Collided s3 at step {}".format(step))
                self.dt_check = self.dt/2.
                # nump = 100

                # If the collision does not happen to be within tolerance
                if abs((self.distance_check - (self.r1 + self.r2))) > self.collision_tol:
                    print("Needed to find the collision position")
                
                    self.dt_change = self.dt_check
                    for a in range(self.max_collision_check):

                        if abs((self.distance_check - (self.r1 + self.r2))) > self.collision_tol:
                            self.precollision = gausss1(startvec = self.simdatalist[step-1],
                                                        params   = self.system_params,
                                                        dynfunc  = dynfunc_h3sim2ballcol,
                                                        dynjac   = dynjac_h3sim2ballcol,
                                                        dt       = self.dt_change,
                                                        tol      = self.tol)
                            
                            self.distance_check = h3dist(rot2hyp(self.precollision[0:3]),rot2hyp(self.precollision[3:6]))

                            # Still overlaping
                            if self.distance_check - (self.r1 + self.r2) < 0.:
                                self.dt_change = self.dt_change - self.dt_check #- dt_check*abs((distck3 - (r1 +r2))/(r1 +r2))/dt
                                print("dt_check minus")
                                print(self.dt_check)
                            # Not overlapping
                            elif self.distance_check - (self.r1 + self.r2) > 0.:
                                self.dt_check = self.dt_check/2
                                self.dt_change = self.dt_change + self.dt_check #+ dt_check*abs((distck3 - (r1 +r2))/(r1 +r2))/dt
                                print("dt_check plus")
                                print(self.dt_check)
                            print('separation')
                            print(self.distance_check - (self.r1 + self.r2))
                            # print(a)
                        else:
                            break
                

                    self.postcollision = h3collision(self.precollision,self.system_params)

                    self.simdatalist[step] = gausss1(startvec   = self.postcollision,
                                                       params   = self.system_params,
                                                       dynfunc  = dynfunc_h3sim2ballcol,
                                                       dynjac   = dynjac_h3sim2ballcol,
                                                       dt       = self.dt - self.dt_change,
                                                       tol      = self.tol)

                # If the collision happens to be within tolerance
                else:
                    print("Collision occurred within tolerance")

                    self.postcollision = h3collision(self.simdatalist[step-1],self.system_params)


                    self.simdatalist[step] = gausss1(startvec   = self.postcollision,
                                                     params     = self.system_params,
                                                     dynfunc    = dynfunc_h3sim2ballcol,
                                                     dynjac     = dynjac_h3sim2ballcol,
                                                     dt         = self.dt,
                                                     tol        = self.tol)
            else:
                # print("Not collided")
                pass


    def set_initial_conditions(self, system_ics):
        """
        Inputs user given initial conditions of system for simulation

        Position and velocity information should be given in terms of the
        parameterization of the ambient space

        Parameters
        ----------
        system_params : array
            the array of parameters describing system consisting of:
                * p1 = initial posiiton of vertex 1
                * p2 = initial posiiton of vertex 2
                * v1 = initial velocity of vertex 1
                * v2 = initial velocity of vertex 2

        """

        self.system_ics = system_ics
        self._have_ics = True

    def clear_data(self):
        """
        Clears any simulation data stored in simulation object

        """

        self.simdatalist = np.zeros((self.t_arr.shape[0],12))
        self._have_run = False

    def run(self):
        """
        Runs simulation once given all necessary information

        Raises
        ----------
        NotImplementedError
            If no initial conditions have been provided
        """

        if self._have_ics:
            # Hyperbolic Space (Sim)
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                self.simdatalist[0] = self.system_ics.copy()

                # Gauss 1-Step Method
                if self.solver_id == "gs1":
                    for step in range(self.t_arr.shape[0] - 1):
                        self.check_for_collision(step)
                        self.simdatalist[step+1] = gausss1(
                            startvec    = self.simdatalist[step],
                            params      = self.system_params,
                            dynfunc     = dynfunc_h3sim2ballcol,
                            dynjac      = dynjac_h3sim2ballcol,
                            dt          = self.dt,
                            tol         = self.tol)
                        
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
            self._have_run = True
        else:
            raise NotImplementedError("Must provide initial conditions via set_initial_conditions() before running simulation")


    def output_data(self):
        """
        Outputs simulation data to file with name:
            {self.ambient_geo}_r_{self.solver_id}_sim_tmax{self.tmax}_dt{self.dt}.npy

        Raises
        ----------
        NotImplementedError
            If simulation has not been run, i.e. no data generated
        """

        if self._have_run:
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                np.save("h3_c_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                np.save("s3_c_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
        else:
            raise NotImplementedError("Must use run() to generate data")

    