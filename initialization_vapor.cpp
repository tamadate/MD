#include "md.hpp"
#include <random>
#include <algorithm>
/*########################################################################################

-----Initialization-----
intialization:
This intialize the calculation. Reading initial positions and setting initial velocities.
reintialization:
This makes connection between thermal relaxation and main diffusion coeficient calculation. It reset position, time, pair list and margine length.

#######################################################################################*/

/////////////////////////////////////////////////////////////////////
/*	
	- Randomly arraying gas molecule around an ion with avoiding the 
	overlapping. The velocity is picked from  the Maxwell-Boltzumann
	distribution
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_vapor(void) {
	double nion=vars->ions.size();
	double vapor[Nof_vapor][3], dx,dy,dz, d,min_gv, min_iv, min_vv, gv, iv, vv;
	gv=10, iv=10, vv=10;	/*	minimum gas-gas, gas-ion distance */
    
    // Set random number generator
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	mt19937 mt(seed());
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);
    
    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	for(int i=0;i<Nof_vapor;i++)	vapor[i][0] = vapor[i][1] = vapor[i][2] = 0;
	int i=0;
	do {
		min_gv=min_iv=min_vv=10000.0;
		vapor[i][0]=r(mt), vapor[i][1]=r(mt), vapor[i][2]=r(mt);
		for(int j=0;j<i;j++) {
			dx = vapor[i][0] - vapor[j][0];
			dy = vapor[i][1] - vapor[j][1];
			dz = vapor[i][2] - vapor[j][2];
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_vv) min_vv=d; // minimum gas-gas distance
		}
		for(auto &a : vars->ions) {
			dx=vapor[i][0]-a.qx;
			dy=vapor[i][1]-a.qy;
			dz=vapor[i][2]-a.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_iv) min_iv=d; // minimum gas-ion distance
		}
        for(auto &a : vars->gases) {
			dx=vapor[i][0]-a.qx;
			dy=vapor[i][1]-a.qy;
			dz=vapor[i][2]-a.qz;
			adjust_periodic(dx, dy, dz);
			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<min_gv) min_gv=d; // minimum gas-ion distance
		}
		if(min_gv>gv && min_iv>iv && min_vv>vv){
			/*double vx=distgas(engine)*1e-5;
			double vy=distgas(engine)*1e-5;
			double vz=distgas(engine)*1e-5;
			if (gastype==1)	vars->add_gases(nion+i+1, 10, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 4.0027);
            if (gastype==2)  vars->add_gases(nion+i+1, 11, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 28.02);
			if (gastype==3)	vars->add_gases(nion+i+1, 12, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 28.02);
			if (gastype==4)	vars->add_gases(nion+i+1, 19, gas[i][0], gas[i][1], gas[i][2], vx, vy, vz, 0.0, 0.0, 0.0, 0.0, 39.948);
			i++;*/
		}
        i++;
	} while(i<Nof_vapor);
    if (gastype==2){
        for(i=0; i<Nof_around_gas; i++){
            vars->add_gases(nion+i+1+Nof_around_gas, 11, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0);
        }
    }
}


