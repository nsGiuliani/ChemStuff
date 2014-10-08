
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
	double temp;
	double dx=0.0;
	double dy=0.0;
	double dz=0.0;
	double box =1000;
	double E = 1.712*(.000000000000000000001);
	double O = 342;
	double Fx;
	double Fy;
	double Fz;
	srand(time(0));
	for(double i =0; i<10; i++) {
		atoms.push_back(Atom(39.948,fmod(rand(),box),fmod(rand(),box),fmod(rand(),box),fmod(rand(),10)-5,fmod(rand(),10)-5,fmod(rand(),10)-5));
	}
	ofstream fout;
	fout.open("test3.xyz");
	for (int i=0; i<10000; i++ ) {
		fout << "10" << endl;
		fout << "iteration " << i << endl;
		for (int k=0; k<10; k++) {
			fout<< "Argon	";
			fout<< atoms[k].getxCoor() <<"	"<<atoms[k].getyCoor() << "	"<< atoms[k].getzCoor()<< endl;
		}
		for (int k=0; k<10; k++) {
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
			for (int l = 0; l<10;l++){
				double min = 1000000000;
				if (l==k) {
					break;
				}
				else{
					temp = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());
					if (temp < min) {
						min = temp;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}
					temp = (atoms[l].getxCoor()-atoms[k].getxCoor()-box)*(atoms[l].getxCoor()-atoms[k].getxCoor()-box)
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (temp < min) {
						min = temp;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor()-box;
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}			
					temp = (atoms[l].getxCoor()-atoms[k].getxCoor()+box)*(atoms[l].getxCoor()-atoms[k].getxCoor()+box)
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (temp < min){
						min = temp;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor()+box;
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}			
					temp = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor()-box)*(atoms[l].getyCoor()-atoms[k].getyCoor()-box)
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (temp < min){
						min = temp;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor()-box;
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}
					temp = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor()+box)*(atoms[l].getyCoor()-atoms[k].getyCoor()+box)
						+ (atoms[l].getzCoor()-atoms[k].getzCoor())*(atoms[l].getzCoor()-atoms[k].getzCoor());					
					if (temp < min){
						min = temp;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor()+box;
						dz =atoms[l].getzCoor()-atoms[k].getzCoor();
					}		
					temp = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor()-box)*(atoms[l].getzCoor()-atoms[k].getzCoor()-box);					
					if (temp < min){
						min = temp;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor()-box;
					}
					temp = (atoms[l].getxCoor()-atoms[k].getxCoor())*(atoms[l].getxCoor()-atoms[k].getxCoor())
						+(atoms[l].getyCoor()-atoms[k].getyCoor())*(atoms[l].getyCoor()-atoms[k].getyCoor())
						+ (atoms[l].getzCoor()-atoms[k].getzCoor()+box)*(atoms[l].getzCoor()-atoms[k].getzCoor()+box);					
					if (temp < min){
						min = temp;
						dx =atoms[l].getxCoor()-atoms[k].getxCoor();
						dy =atoms[l].getyCoor()-atoms[k].getyCoor();
						dz =atoms[l].getzCoor()-atoms[k].getzCoor()+box;
					}
					// input code here for updating the velocity
					double D = dx*dx+dy*dy+dz*dz;
					double D_7= (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)
						*(dx*dx+dy*dy+dz*dz);
					double D_4 = (dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz)*(dx*dx+dy*dy+dz*dz);
					double o_6= O*O*O*O*O*O;
					Fx = (24*E*o_6*dx*(2*o_6/(D_7)+1/(D_4)));
					Fy = (24*E*o_6*dy*(2*o_6/(D_7)+1/(D_4)));					
					Fz = (24*E*o_6*dz*(2*o_6/(D_7)+1/(D_4)));
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
