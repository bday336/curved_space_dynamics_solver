# // a class for integrating equations of motion called "integrator" and one specific implementation, Runge Kutta
# //derive is a function taking a state to state (now storing velocity and acceleration instead of position and velocity)
# //items fed into RungeKutta need to have the following methods available:
# //.add, .multiplyScalar, .clone

# //implementing the Rk4 Scheme for arbitrary classes that have clone add and multiplyScalar
# //will use this possibly on individual states, or on entire DataLists!
class RungeKutta:

    def __init__(self, derive, ep):
        self.derive=derive
        self.ep=ep

    # //step forwards one timestep
    def step(self, state):

        # k1,k2,k3,k4;
        # temp;

        # //get the derivative
        k1 = self.derive(state)
        k1.multiplyScalar(self.ep)

        # //get k2
        temp=state.clone().add(k1.clone().multiplyScalar(0.5))
        k2=self.derive(temp)
        k2.multiplyScalar(self.ep)

        # //get k3
        temp=state.clone().add(k2.clone().multiplyScalar(0.5))
        k3=self.derive(temp)
        k3.multiplyScalar(self.ep)

        # //get k4
        temp=state.clone().add(k3.multiplyScalar(1.))
        k4=self.derive(temp)
        k4.multiplyScalar(self.ep)

        # //add up results:
        total = k1 #//scale factor 1
        total.add(k2.multiplyScalar(2))
        total.add(k3.multiplyScalar(2))
        total.add(k4) #//scale factor 1
        total.multiplyScalar(1/6)

        # //move ahead one step
        return state.clone().updateBy(total)