
#include<iostream>
#include<cmath>
#include<cstdlib>//needed for rand
#include<ctime> //needed for srand(time(0))
#include<fstream> //needed for files
#include<vector>
#include "atom.cpp"

using namespace std;

int main() {

	vector<Atom> atoms;
	double x;
	double y;
	double z;
	double dist;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	double box =727;
	double E = 1.712*(.000000000000000000001);
	double sig = 3.4*(.0000000001);
	double Fx;
	double Fy;
	double Fz;
	int numOfAtoms=10;
	double timeStep = 1.0;
	int iterations = 10000;
	srand(time(0));
	//set initial positions and velocities for atoms
	for(double i =0; i<numOfAtoms; i++) {
		atoms.push_back(Atom(39.948,fmod(rand(),box),fmod(rand(),box),fmod(rand(),box),fmod(rand(),50)-25,fmod(rand(),50)-25,fmod(rand(),50)-25));
	}
	ofstream fout;
	fout.open("10Atom.xyz");
	//print to .xyz file
	for (int i=0; i<iterations; i++ ) {
		fout << numOfAtoms << endl;
		fout << "iteration " << i << endl;
		for (int k=0; k<numOfAtoms; k++) {
			fout<< "Argon	";
			fout<< atoms[k].getxCoor() <<"	"<<atoms[k].getyCoor() << "	"<< atoms[k].getzCoor()<< endl;
		}
		for (int k=0; k<numOfAtoms; k++) { //k is the atom which is being updated
			//update the coordinates and make sure they are in the box
			x= atoms[k].getxCoor()+atoms[k].getxVel();
			y= atoms[k].getyCoor() + atoms[k].getyVel();
			z=  atoms[k].getzCoor() + atoms[k].getzVel();
			if(x >= box){
				x =fmod(x, box);
			}
			if(x <0) {
				x =fmod(x, box) + box;
			}
			if(y >= box){
				y = fmod(y, box);
			}
			if(y <0) {
				y =fmod(y, box) +box;
			}
			if(z >= box){
				z =fmod(z, box);
			}
			if(z <0) {
				z =fmod(z, box) + box;
			}
			atoms[k].set_Coor(x,y,z);
			//find the minimum distance to other atoms and update vel
			for (int l = 0; l<numOfAtoms;l++){ //l is the atom whose presence is applying force to k
				double min = 1000000000;
				if(l != k){ //make sure they are not the same atom
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());
					if (dist < min) {
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor()-box)*(atoms[l].getxCoor()-atoms[k].getxCoor()-box)
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min) {
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor()-box;
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}			
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor()+box)*(atoms[l].getxCoor()-atoms[k].getxCoor()+box)
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor()+box;
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}			
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor()-box)*(atoms[l].getyCoor()-atoms[k].getyCoor()-box)
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor()-box;
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor()+box)*(atoms[l].getyCoor()-atoms[k].getyCoor()+box)
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor()+box;
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}		
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor()-box)*(atoms[l].getzCoor()-atoms[k].getzCoor()-box);					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor()-box;
					}
					dist = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor()+box)*(atoms[l].getzCoor()-atoms[k].getzCoor()+box);					
					if (dist < min){
						min = dist;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor()+box;
					}
					// update the velocity 
					double D = dx*dx+dy*dy+dz*dz;
					double D_7= (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)
						*(dx*dx+dy*dy+dz*dz);
					double D_4 = (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz);
					double sig_6= sig*sig*sig*sig*sig*sig;
					Fx = (24*E*sig_6*dx*(2*sig_6/(D_7)+1/(D_4)))*timeStep;
					Fy = (24*E*sig_6*dy*(2*sig_6/(D_7)+1/(D_4)))*timeStep;					
					Fz = (24*E*sig_6*dz*(2*sig_6/(D_7)+1/(D_4)))*timeStep;
					x= atoms[k].getxVel()+Fx/atoms[k].getMass();
					y= atoms[k].getyVel() + Fy/atoms[k].getMass();
					z=  atoms[k].getzVel() + Fz/atoms[k].getMass();
					atoms[k].set_Vel(x,y,z);
									
				}		
																	
			}

		}
		
	}

	return 0;
}
