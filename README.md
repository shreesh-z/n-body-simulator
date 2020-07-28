# n-body-simulator
An N-Body Simulator I'm building using the B-H Tree algorithm.
Right now, it is able to build a B-H quad tree from an array of body objects, and find the forces acting on them by traversing through it.
It is also able to evolve the system in time using an RK4 integrator.
Planned:
* Real time rendering of the system with openGL

Once the C++ build is finished, I'm also planning on implementing the simulator with python (numpy+pygame) and/or Julia.
