# Dynamics Solver for Extended Body Systems in Generic Curved Space

Source code for Generic Curved Space Dynamics Solver developed as part of my PhD thesis work in collaboration with Steve Trettel and Sabetta Matsumoto.

Current version of simulation scripts found here is designed with the constant curvature spaces space (3-dimensional spherical S^3, euclidean E^3, and hyperbolic H^3 spaces in mind). Also included are the test simulations for a non-homogeneous, isotropic space with an analytic distance function used to test "curvature detector" system. Simulations involving harmonic coupling potential as well as rigid rod constraints between discrete elements of extended body system can be implemented in constant curvature spaces for two main test systems: rod body and triangle body. Collisions are facilitated via hard sphere elastic collisions between the mass vertices making up the discretized extended body system.

Progress on this code is ongoing to provide a production ready dynamics solver. This work is also being developed in parallel with a academic publication where we discuss the theoretical models we use for generic curved spaces as well as present the results of the analysis for test systems in the constant curvature geometries.
