# n-body-simulator
An N-Body Simulator I'm building using the B-H Tree algorithm.

It builds a B-H quad tree from an array of body objects, and finds the forces acting on them by traversing through it.
Then it evolves the system in time using an RK4 integrator.

The interface uses openGL for rendering. For this reason, a mingw32 g++ compiler has to be used to build the project.
The .dll binary files are needed for the .exe file to run properly.

The barnes_hut.cpp file hosues the structs, data structures and the algorithm needed to represent the system and its ODE.
The nbody-simulator.cpp file houses the rendering. The two files interact using the barnes_hut.h header file.

Now that the C++ build is finished, I'm also planning on implementing the simulator with python (numpy+pygame) and/or Julia.
