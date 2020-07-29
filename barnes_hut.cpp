#include <stdio.h>
#include <random>
#include <cmath>

#include "barnes_hut.h"

/*
	------------------------GLOBAL VARIABLES---------------------------
*/

//number of bodies in the space
int N = 10;
//Gravitation constant
double G = 0.0000000005;
//theta criterion (inverse)
double THETA_CRIT_INV = 2.0;
//Cap on the proportionality coefficient of the grav force
double GRAV_COEFF_CAP = 10000.0;
//half dimensions of the square space
double SQUARE_DIMS = 1000.0;
//bounds of the space in cartesian coordinates
double LEFT_POS = -SQUARE_DIMS;
double DOWN_POS = -SQUARE_DIMS;
double HALF_SIDE_LENGTH = SQUARE_DIMS;

/*
	---------------------------BODY STRUCT-----------------------------
	
	Attributes:
		ID     : 0 for interior node bodies (COM bodies) and positive int for real bodies
		mass   : mass of the body
		coords : x,y position and vx,vy velocity of the body
		force  : x,y gravitational force acting on the body
		
	Constructors:
		Body(ID, m_, posX, posY, velX, velY)
			for a general body
		Body(ID, m_, posX, posY)
			for a body starting from rest
		Body()
			for an interior node in the B-H Tree
	
	Member Functions:
		updateCOM(Body*):
			For bodies of interior nodes with ID = 0. When a new body is added to the subtree
			of the interior node, this function is called on the COM body of the node.
		initForce():
			Initializes the force acting on the body to 0.
*/

struct Body{
	//ID is used while calculating forces
	//If the ID of the body is the same as of the body of
	//the external node, force is zero
	int ID;
	double mass;
	
	//ordering is posX, posY, velX, velY;
	double coords[4];
	
	//Force acting on the body
	double force[2];
	
	//For a general body
	Body(int ID_, double m_, double posX_, double posY_, double velX_, double velY_){
		ID = ID_;
		mass = m_;
		coords[0] = posX_;
		coords[1] = posY_;
		coords[2] = velX_;
		coords[3] = velY_;
		force[0] = force[1] = 0.0;
	}
	
	//For a body starting from rest
	Body(int ID_, double m_, double posX_, double posY_) 
	: Body(ID_, m_, posX_, posY_, 0.0, 0.0){}
	
	//Empty body
	Body()
	: Body(0, 0.0, 0.0, 0.0){}
	
	void updateCOM(Body *otherBod){
		//Updates COM with other body's attributes
		//Other body's attributes remain unchanged
		double totMass = mass + otherBod->mass;
		for( int i = 0; i < 4; i++ ){
			coords[i] = ( coords[i]*mass + otherBod->coords[i]*otherBod->mass )/totMass;
		}
		mass = totMass;
	}
	
	//initializes force on body to zero
	void initForce(){
		force[0] = force[1] = 0.0;
	}
};

//The array of bodies that are being simulated
Body **Bodies;


/*
	-------------------------BARNES-HUT TREE NODE STRUCT--------------------------------
	
	Attributes:
		*body      : pointer to the body contained in the node
		left, down : coordinates of the bottom-left corner of the square space represented by the node
		s          : half side length of the square represented by the node
		*quadNodes : pointer array that houses the 4 sub trees of the node
		
	Constructors:
		Node(*body, left_, down_, side)
			for a general node with a general body contained in it
		Node(left_, down_, side)
			for an internal node that starts with having at rest with no mass
*/

struct Node{
	//Either the body inside the node or the COM body
	Body *body;
	//The left and down coords of the quad tree cell
	double left, down;
	//The HALF of side length of the cell (NOTE : HALF OF SIDE LENGTH)
	double s;
	
	//The four nodes inside the given node
	//Since four nodes, a quad tree is formed
	//node ordering: NE, SE, NW, SW
	Node* quadNodes[4];
	
	//Constructor for a general node
	Node(Body *body_, double left_, double down_, double side){
		body = body_;
		s = side;
		left = left_;
		down = down_;
		for( int i = 0; i < 4; i++){
			quadNodes[i] = NULL;
		}
	}
	//Constructor for an internal node, having ID 0
	Node(double left_, double down_, double side)
	: Node(new Body(0, 0.0, 0.0, 0.0), left_, down_, side){}
};

/*
	---------------------------------FUNCTIONS----------------------------------
	
	insert:
		Inserts a body into the BH Tree recursively
	printTree:
		prints the BH Tree recursively (for debugging only)
	findForce:
		finds the force on a body due to all the other bodies inside the BHTree
	deleteTree:
		deletes the entire BHTree to start fresh in the next iteration
	RK4:
		RK4 integrator of the whole system (4xN dimension ODE solver)
	init:
		initializes the Bodies array
	evolve:
		evolves the bodies in time
	main:
		main function for running the program on the console
	
*/

//Inserts the body inBody into the B-H quad tree recursively
void insert(Node *root, Body *inBody){
	if( root->body->ID == 0 ){
		//If the root node is an interior node
		
		//The index of one of four nodes in which the newly added body is located
		//0: NE, 1:SE, 2:NW, 3:SW
		int nodeListInd = 0;
		
		//To store offsets away from the current node's (left, down)
		//where a new interior node will be created
		int directionX = 0, directionY = 0;
		
		//Updating COM of combined system
		root->body->updateCOM(inBody);
		
		//now figure out which direction quad node to place body in
		if( inBody->coords[0] >= (root->left+root->s) ){
			//If body is in East
			if( inBody->coords[1] >= (root->down+root->s) ){
				//If body is in NE
				nodeListInd = 0;
				//Take steps in both directions
				directionX = 1;
				directionY = 1;
			}else{
				//If body is in SE
				nodeListInd = 1;
				//Take step only in x direction
				directionX = 1;
			}
		}else{
			//If body is in West
			if( inBody->coords[1] >= (root->down+root->s) ){ 
				//If body is in NW
				nodeListInd = 2;
				//Take step only in y direction
				directionY = 1;
			}else{
				//If body is in SW
				nodeListInd = 3;
				//No steps to take
			}
		}
		if( root->quadNodes[nodeListInd] == NULL ){
			//Make new node if null, then exit
			root->quadNodes[nodeListInd] = new Node(inBody, 
													root->left + directionX * root->s, 
													root->down + directionY * root->s, 
													root->s/2);
			return;
		}else{
			//If a node exists there already, add recursively
			insert(root->quadNodes[nodeListInd], inBody);
			return;
		}
	}else{
		//If root node is an exterior node
		
		//Making a new interior COM body to replace the exterior node
		Body *newBod = new Body();
		
		//storing the address of the old exterior node's body
		Body *oldBod = root->body;
		
		//Setting the body of the node to the interior node-body
		root->body = newBod;
		
		//Recursively add the two bodies into the root node again
		insert(root, oldBod);
		insert(root, inBody);
		return;
	}
}

//print the B-H tree in post order
void printTree(Node *root){
	for(int i = 0; i < 4; i++){
		if(root->quadNodes[i] != NULL){
			printTree(root->quadNodes[i]);
		}
	}
	//std::cout << root->body->ID << ' ' << root->body->mass << ' ' << root->s << "\n";
	printf("%d %lf %lf\n", root->body->ID, root->body->mass, root->s);
}

//find net force on a body due to the bodies present inside the B-H tree
//acts recursively
void findForce(Body *body, Node *root){
	
	//If the same body is encountered in the tree, don't add anything to the force
	if( body->ID == root->body->ID )
		return;
	
	//X and Y distances between the two bodies
	double diffX = root->body->coords[0] - body->coords[0];
	double diffY = root->body->coords[1] - body->coords[1];
	
	//Net distance between two bodies
	double dist = std::hypot(diffX, diffY);
	
	if( dist > THETA_CRIT_INV*root->s*2 || root->body->ID != 0){
		//If COM body satisfies theta criterion
		//OR if an exterior node has been reached
		double coeff = (G * body->mass * root->body->mass)/std::pow(dist,3);
		
		//capping coeff value so infinities are not produced
		coeff = coeff > GRAV_COEFF_CAP ? GRAV_COEFF_CAP : coeff;
		
		//calculating force acc. to Newton's Law of Gravitation
		body->force[0] += coeff*diffX;
		body->force[1] += coeff*diffY;
		return;
	}
	else{
		//COM body doesn't satisfy theta criterion, go deeper into the tree
		//iterating over all the 4 nodes inside current node
		for( int i = 0; i < 4; i++ ){
			if( root->quadNodes[i] != NULL ){
				//recursively finding the force on current body by bodies inside current node
				findForce(body, root->quadNodes[i]);
			}
		}
		return;
	}
}

//to delete tree once all forces have been found
void deleteTree(Node *node){
	if( node == NULL )
		return;
	for( int i = 0; i < 4; i++ )
		deleteTree(node->quadNodes[i]);
	
	delete node;
}

//RK4 integrator
void RK4(Body **Bodies, double timestep){
	
	//Coords of bodies at time t_init
	double *OrigCoords = new double[N*4]; //Nx4 matrix
	for( int i = 0; i < N; i++){
		for( int j = 0; j < 4; j++){
			OrigCoords[4*i + j] = Bodies[i]->coords[j];
		}
	}
	
	//stores the 4 sets of RK4 function values
	double *functs = new double[N*4*4];
	
	for( int rk = 0; rk < 4; rk++ ){
		
		//The root node of the BH tree, used to find the forces
		Node *root = new Node(LEFT_POS, DOWN_POS, HALF_SIDE_LENGTH);
		
		//First add all bodies to the BH Tree
		for( int i = 0; i < N; i++ ){
			//Initialize forces with 0
			Bodies[i]->initForce();
			insert(root, Bodies[i]);
		}
		
		//Then find force acting on each body
		for( int i = 0; i < N; i++ ){
			findForce(Bodies[i], root);
		}
		
		//Now delete the tree for a fresh tree in the next RK4 iteration
		deleteTree(root);
		
		//storing function values found at this RK step (ie. velocity and force)
		for( int i = 0; i < N; i++ ){
			for( int j = 0; j < 2; j++ ){
				//storing velocity values as functs for displacement
				functs[4*i + j + 4*N*rk] = Bodies[i]->coords[j+2];
				//storing force values as functs for velocity
				functs[4*i + (j+2) + 4*N*rk] = Bodies[i]->force[j];
			}
		}
		
		//Now forces on all bodies have been found. RK4 Step is to be taken.
		
		if( rk < 2 ){
			for( int i = 0; i < N; i++ ){
				for( int j = 0; j < 2; j++ ){
					//First x and y coords
					Bodies[i]->coords[j] = OrigCoords[4*i + j] + (timestep/2)*Bodies[i]->coords[j+2];
					//Now vx and vy
					Bodies[i]->coords[j+2] = OrigCoords[4*i + j+2] + (timestep/2)*Bodies[i]->force[j];
				}
			}
		}else if( rk == 2 ){
			for( int i = 0; i < N; i++ ){
				for( int j = 0; j < 2; j++ ){
					//First x and y coords
					Bodies[i]->coords[j] = OrigCoords[4*i + j] + timestep*Bodies[i]->coords[j+2];
					//Now vx and vy
					Bodies[i]->coords[j+2] = OrigCoords[4*i + j+2] + timestep*Bodies[i]->force[j];
				}
			}
		}else{
			//final RK4 step
			for( int i = 0; i < N; i++ ){
				for( int j = 0; j < 4; j++ ){
					//taking full RK4 step with the 4xNx4 matrix
					Bodies[i]->coords[j] = OrigCoords[4*i + j] + (timestep/6)*(
																		functs[4*i + j] +
																		2*( functs[4*i + j + 4*N] + functs[4*i + j + 8*N] ) +
																		functs[4*i + j + 12*N]
																		);
				}
				//Bounding the x and y coords in the square
				if( Bodies[i]->coords[0] < -HALF_SIDE_LENGTH )
					Bodies[i]->coords[0] += 2*HALF_SIDE_LENGTH;
				else if( Bodies[i]->coords[0] > HALF_SIDE_LENGTH )
					Bodies[i]->coords[0] -= 2*HALF_SIDE_LENGTH;
				
				if( Bodies[i]->coords[1] < -HALF_SIDE_LENGTH )
					Bodies[i]->coords[1] += 2*HALF_SIDE_LENGTH;
				else if( Bodies[i]->coords[1] > HALF_SIDE_LENGTH )
					Bodies[i]->coords[1] -= 2*HALF_SIDE_LENGTH;
			}
		}
	}
}

void initBodies(){
	//pRNG functions from stdlib
	std::mt19937 gen(1);		//constant seed during development
	std::uniform_real_distribution<> disMass(1.0, 2.0), disPos(LEFT_POS/2,-LEFT_POS/2);		
	
	//filling the bodies array with random masses and initial positions
	for(int i = 0; i < N; i++){
		//Bodies are ID'd starting from 1
		Bodies[i] = new Body(i+1, disMass(gen), disPos(gen), disPos(gen));
	}
}

void evolveBodies(){
	double maxTime = 5;
	double timestep = 1;
	int count = (int)(maxTime/timestep);
	
	for( int t = 0; t < count; t++ ){
		RK4(Bodies, timestep);
		for( int i = 0; i < N; i++){
			//std::cout << Bodies[i]->coords[0] << ' ' << Bodies[i]->coords[1] << ' ' << Bodies[i]->coords[2] << ' ' << Bodies[i]->coords[2] << "\n";
			printf("%f\n",Bodies[i]->coords[0]);
		}
		printf("\n");
		//std::cout << "\n";
	}
}

void evolveBodiesStep(double timestep){
	RK4(Bodies, timestep);
}

void setSimulationParams(int N_, double G_, double theta_, double coeffCap_, double squareDims_){
	N = N_;
	G = G_;
	THETA_CRIT_INV = theta_;
	GRAV_COEFF_CAP = coeffCap_;
	SQUARE_DIMS = squareDims_;
	Bodies = new Body*[N];
}

double* getBodyPoints(){
	//Nx2 matrix with x,y coordinates of the bodies
	double *points = new double[N*2];
	for( int i = 0; i < N; i++ ){
		points[2*i] = Bodies[i]->coords[0];
		points[2*i + 1] = Bodies[i]->coords[1];
	}
	return points;
}

/*int main(){
	//array holding all the bodies
	//Body *Bodies[N];
	initBodies();
	evolveBodies();
	return 0;
}*/
