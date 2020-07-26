#include <iostream>
#include <random>
#include <cmath>

//for value of gravitation constant
#define G 0.0000000000005

struct Body{
	/*ID is used while calculating forces
	If the ID of the body is the same as of the body of
	the external node, force is zero*/
	int ID;
	double mass;
	double coords[4];
	//ordering is posX, posY, velX, velY;

	//For a general body
	Body(int ID_, double m_, double posX_, double posY_, double velX_, double velY_){
		ID = ID_;
		mass = m_;
		coords[0] = posX_;
		coords[1] = posY_;
		coords[2] = velX_;
		coords[3] = velY_;
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
			coords[i] = (coords[i]*mass + otherBod->coords[i]*otherBod->mass )/totMass;
		}
		mass = totMass;
	}
};

struct Node{
	//Either the body inside the node or the COM body
	Body *body;
	//The left and down coords of the quad tree cell
	double left, down;
	//The half of side length of the cell
	double s;
	
	//The four nodes inside the given node
	//Since four nodes, a quad tree is formed
	Node* quadNodes[4];
	//node ordering: SE, NE, SW, NW
	
	//Constructor for a general body
	Node(Body *bod_, double left_, double down_, double side){
		body = bod_;
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

//Inserts the body inBody into the B-H quad tree recursively
void insert(Node *root, Body *inBody){
	if( root->body->ID == 0 ){
		//If the root node is an interior node
		
		//The index of list of quad nodes
		int nodeListInd = 0;
		//To store offsets if a new node has to be created
		int directionX = 0, directionY = 0;
		
		//First update COM of combined system
		root->body->updateCOM(inBody);
		
		//now figure out which quad node to place body in
		if( inBody->coords[0] >= (root->left+root->s) ){
			//If body is in East
			if( inBody->coords[1] >= (root->down+root->s) ){
				//If body is in NE
				nodeListInd = 0;
				//Take steps in both direction
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
		
		//Making a new interior COM node to replace the exterior node
		Body *newBod = new Body();
		Body *rootBod = root->body;
		
		//Adding the ext node body and new inBody to the new interior node
		newBod->updateCOM(rootBod);
		newBod->updateCOM(inBody);
		
		//Setting the body of the node to an interior node
		root->body = newBod;
		
		//Recursively add the two bodies into the root node again
		insert(root, rootBod);
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
	std::cout << root->body->ID << ' ' << root->body->mass << ' ' << root->s << "\n";
}

//Used to find steps during RK4, unused now
/*
double* step(double *functs, double *coords, double timestep){
	double newVect[4];
	for( int i = 0; i < 4; i++ ){
		newVect[i] = coords[i] + functs[i]*timestep;
	}
	return newVect;
}*/

//returns the net force on body due to the bodies present inside the B-H tree
//acts recursively
//call this function with an initial force {0.0, 0.0}
double* findForce(double* force, Body *body, Node *root, double thetaCritInv){
	
	//If the same body is encountered in the tree, don't add anything to the force
	if( body->ID == root->body->ID )
		return force;
	
	//X and Y distances between the two bodies
	double diffX = body->coords[0] - root->body->coords[0], diffY = body->coords[1] - root->body->coords[1];
	
	//Net distance between two bodies
	double dist = std::hypot(diffX, diffY);
	
	if( dist > thetaCritInv*root->s*2 ){
		//If COM body satisfies theta criterion
		double coeff = (G * body->mass * root->body->mass)/std::pow(dist,3);
		force[0] += coeff*diffX;
		force[1] += coeff*diffY;
		return force;
	}
	else{
		//iterating over all the 4 nodes inside current node
		for( int i = 0; i < 4; i++ ){
			if( root->quadNodes[i] != NULL ){
				//finding the force on current body by bodies inside current node recursively
				double *newForce = new double[2];
				newForce = findForce(force, body, root->quadNodes[i], thetaCritInv);
				force[0] += newForce[0];
				force[1] += newForce[1];
			}
		}
		return force;
	}
}
 

int main(){
	//Total number of bodies
	int N = 10;
	
	//half dimensions of the square space
	double dims = 1000.0;
	double leftPos = -dims, downPos = -dims, halfSide = dims;
	int thetaCritInv = 2;	//In this case, theta criterion is 0.5, ie d > 2*sidelength
	
	//Root node of the barnes hut quad tree
	Node* root = new Node(leftPos, downPos, halfSide);
	
	//array holding all the bodies
	Body* Bodies[N];
	
	//pRNG functions from stdlib
	std::mt19937 gen(1);		//constant seed during development
	std::uniform_real_distribution<> disMass(1.0, 2.0), disPos(leftPos,-leftPos);		
	
	//filling the bodies array with random masses and initial positions
	for(int i = 0; i < N; i++){
		//Bodies are ID'd starting from 1
		Bodies[i] = new Body(i+1, disMass(gen), disPos(gen), disPos(gen));
	}
	
	//Adding bodies to the BH tree
	for (int i = 0; i < N; i++){
		insert(root, Bodies[i]);
	}
	
	//std::cout << "Tree traversal:\n";
	//printTree(root);
	//std::cout << "\n";
	
	std::cout << "Forces on bodies:\n";
	for( int i = 0; i < N; i++ ){
		double *force = new double[2];
		force = findForce(force, Bodies[i], root, thetaCritInv);
		std::cout << force[0] << ' ' << force[1] << "\n";
	}
	
	return 0;
}