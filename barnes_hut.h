#ifndef BARNES_H
#define BARNES_H
	void initBodies();
	void evolveBodies();
	void evolveBodiesStep(double timestep);
	void setSimulationParams(int N_, double G_, double theta_, double coeffCap_, double squareDims_);
	double* getBodyPoints();
#endif