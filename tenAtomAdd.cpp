
#include<iostream>
#include<cmath>
#include<cstdlib>//needed for rand
#include<ctime> //needed for srand(time(0))
#include<fstream> //needed for files
#include<vector>
#include "atom.cpp"

using namespace std;

void simulation(int numAtoms, double temperature) {

	vector<Atom> atoms;
	double x;
	double y;
	double z;
	double dist;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;

	double box =0;

	double E = 1.712*(.000000000000000000001);
	double sig = 3.4*(.0000000001);
	double Fx;
	double Fy;
	double Fz;
	int numOfAtoms=10;
	int iterations = 50000;
	double timeStep = .00000000000001;
	double ang = 0.0000000001;
	double ePot = 0;
	double eKin = 0;
	double eTot1;
	double eTot2;
	
	//calculate the volume of the box using PV = nRT (P = 1)
	double n = numAtoms/(6.022*pow(10, 23)); //calculate the mols of atoms
	double R = .08206; //gas constant
	box = n*R*temperature; //calculate the volume of the box in liters
	box *= .001; //convert the volume of the box from liters to meters cubed
	box = pow(box, 1/3.); //take the cube root to find the length of one side of the box
	box /= ang; //convert to angstroms

	srand(time(0));

	//make a box
	ofstream fout;
	fout.open("box.xyz");
	fout << 12*box << endl;
	fout << "box " << endl;
	for (int b =0; b<box; b++) {
		
		fout<<"box   "<< b <<"	"<< 0 << "	"<< 0 << endl;
		fout<<"box   "<< 0 <<"	"<< b << "	"<< 0 << endl;
		fout<<"box   "<< 0 <<"	"<< 0 << "	"<< b << endl;
		fout<<"box   "<< box <<"	"<< b << "	"<< 0 << endl;
		fout<<"box   "<< box <<"	"<< 0 << "	"<< b << endl;
		fout<<"box   "<< b <<"	"<< box << "	"<< 0 << endl;
		fout<<"box   "<< 0 <<"	"<< box << "	"<< b << endl;
		fout<<"box   "<< b <<"	"<< 0 << "	"<< box << endl;
		fout<<"box   "<< 0 <<"	"<< b << "	"<< box << endl;
		fout<<"box   "<< box <<"	"<< box << "	"<< b << endl;
		fout<<"box   "<< b <<"	"<< box << "	"<< box << endl;
		fout<<"box   "<< box <<"	"<< b << "	"<< box << endl;
	}
	fout.close();

	//set initial positions and velocities for atoms
	for(double i =0; i<numOfAtoms; i++) {
		atoms.push_back(Atom(39.948*1.66053892 * pow(10,-27),fmod(rand(),box),fmod(rand(),box),fmod(rand(),box),fmod(rand(),500)-250,fmod(rand(),500)-250,fmod(rand(),500)-250));
	}
	//ofstream fout;
	fout.open("10Atom.xyz");
	//print to .xyz file
	for (int i=0; i<iterations; i++ ) {
		ePot = 0;
		eKin = 0;
		fout << numOfAtoms << endl;
		fout << "iteration " << i << endl;
		for (int k=0; k<numOfAtoms; k++) {
			fout<< "Argon	";
			fout<< atoms[k].getxCoor() <<"	"<<atoms[k].getyCoor() << "	"<< atoms[k].getzCoor()<< endl;
		}
		for (int k=0; k<numOfAtoms; k++) { //k is the atom which is being updated
			//update the coordinates and make sure they are in the box
			x= atoms[k].getxCoor() + atoms[k].getxVel()*timeStep/ang;
			y= atoms[k].getyCoor() + atoms[k].getyVel()*timeStep/ang;
			z= atoms[k].getzCoor() + atoms[k].getzVel()*timeStep/ang;
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
					dx= dx*ang;
					dy=dy*ang;
					dz=dz*ang;
					// update the velocity 
					double D = dx*dx+dy*dy+dz*dz;
					double D_7= (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)
						*(dx*dx+dy*dy+dz*dz);
					double D_4 = (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz);
					double sig_6= sig*sig*sig*sig*sig*sig;
					Fx = -1*(24*E*sig_6*dx*((2*sig_6/(D_7))-(1/(D_4))));
					Fy = -1*(24*E*sig_6*dy*((2*sig_6/(D_7))-(1/(D_4))));					
					Fz = -1*(24*E*sig_6*dz*((2*sig_6/(D_7))-(1/(D_4))));
					x= atoms[k].getxVel() + (Fx/atoms[k].getMass())*timeStep;
					y= atoms[k].getyVel() + (Fy/atoms[k].getMass())*timeStep;
					z= atoms[k].getzVel() + (Fz/atoms[k].getMass())*timeStep;
					atoms[k].set_Vel(x,y,z);
					ePot += 4*E*(pow((sig/sqrt(dist)),12)-pow((sig/sqrt(dist)),6));
				}		
																	
			}
			eKin += .5*atoms[k].getMass()*(sqrt(atoms[k].getxVel()*atoms[k].getxVel()+atoms[k].getyVel()*atoms[k].getyVel()+atoms[k].getzVel()*atoms[k].getzVel()));

		}
		eTot1 = eKin + ePot;
		eTot2 = (eKin + ePot)/1000*6.02*pow(10,23);
		if (i%100 ==0) {
			printf("%d  %4.16g        %4.16g \n",i, eTot1, eTot2);
		}
	}
}

int main() {
	
	double numberOfAtoms;
	double temp;

	cout << "Enter the number of atoms: ";
	cin >> numberOfAtoms;
	cout << "Enter the temperature: ";
	cin >> temp;

	simulation(numberOfAtoms, temp);
	return 0;
}
