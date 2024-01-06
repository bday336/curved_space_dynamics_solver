import numpy as np
from src.integrator_files.integrator_bank import h3rads2tridae, h3rads3tridae, s3rads2tridae, s3rads3tridae
from src.test_system_simulations.test_system_bank import dynfunc_h3simtriangle, dynjac_h3simtriangle, dynfunc_s3simtriangle, dynjac_s3simtriangle


# Simulation class setup

class RigidTriangleSimulation:
    """
    A class used to perform rod body simulation with rigid rod constraints

    ...

    Attributes
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
            * m3 = mass of vertex 3
            * x1 = fixed rigid constraint length (side 12) 
            * x2 = fixed rigid constraint length (side 13) 
            * x3 = fixed rigid constraint length (side 23) 

    dt : float
        the simuation time step size
    tmax : float
        the total simulation time

    solver_id : str
        the solver to be used to evaluate dynamics
        - currently supporting:
            * rs2 = 2-step Radau collocation (RadauIIA based on the work of Schweizer and Li 2015)
            * rs3 = 3-step Radau collocation (RadauIIA based on the work of Schweizer and Li 2015)

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
            {self.ambient_geo}_t_rig_{self.solver_id}_sim_tmax{self.tmax}_dt{self.dt}.npy
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
                * m3 = mass of vertex 3
                * x1 = fixed rigid constraint length (side 12) 
                * x2 = fixed rigid constraint length (side 13) 
                * x3 = fixed rigid constraint length (side 23) 

        dt : float
            the simuation time step size
        tmax : float
            the total simulation time

        solver_id : str
            the solver to be used to evaluate dynamics
            - currently supporting:
                * rs2 = 2-step Radau collocation (RadauIIA based on the work of Schweizer and Li 2015)
                * rs3 = 3-step Radau collocation (RadauIIA based on the work of Schweizer and Li 2015)

        """
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

        self._have_ics = False
        self._have_run = False

    def set_initial_conditions(self, system_ics):
        """
        Inputs user given initial conditions of system for simulation

        Position and velocity information should be given in terms of the
        parameterization of the ambient space

        Parameters
        ----------
        system_params : array
            the array of parameters describing system consisting of:
                * p1   = initial posiiton of vertex 1
                * p2   = initial posiiton of vertex 2
                * p3   = initial posiiton of vertex 3
                * v1   = initial velocity of vertex 1
                * v2   = initial velocity of vertex 2
                * v3   = initial velocity of vertex 3
                * lam1 = initial value of Lagrange multiplier for side 12
                * lam2 = initial value of Lagrange multiplier for side 13
                * lam3 = initial value of Lagrange multiplier for side 23

        """
        self.system_ics = system_ics
        self._have_ics = True

    def clear_data(self):
        """
        Clears any simulation data stored in simulation object

        """
        self.simdatalist = np.zeros((self.t_arr.shape[0],18+3))
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
            self._have_run = True
        else:
            raise NotImplementedError("Must provide initial conditions via set_initial_conditions() before running simulation")

    def output_data(self):
        """
        Outputs simulation data to file with name:
            {self.ambient_geo}_t_rig_{self.solver_id}_sim_tmax{self.tmax}_dt{self.dt}.npy

        Raises
        ----------
        NotImplementedError
            If simulation has not been run, i.e. no data generated
        """
        if self._have_run:
            if self.ambient_geo == "h3" or self.ambient_geo == "H3":
                np.save("h3_t_rig_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
            if self.ambient_geo == "s3" or self.ambient_geo == "S3":
                np.save("s3_t_rig_{}_sim_tmax{}_dt{}".format(self.solver_id, str(self.tmax), str(self.dt)), self.simdatalist)
        else:
            raise NotImplementedError("Must use run() to generate data")

    