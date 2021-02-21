#include "md.hpp"
#include <random>
#include <algorithm>

//////////////////////////////////////////////////////////////////////
/*	
	- Reset time (time -> 0).	
	- Re-initializating the ion and gas positions and velocities. 
	Position -> 0 and delta, Traslational velocity -> 0
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
*/
/////////////////////////////////////////////////////////////////////
void
MD::reinitialization(void) {
	vars->time=0;

	for (auto &a : vars->ions) a.qx-=ion_r[0], a.qy-=ion_r[1], a.qz-=ion_r[2];
	for (auto &a : vars->gases) a.qx-=ion_r[0], a.qy-=ion_r[1], a.qz-=ion_r[2];

	analysis_ion();
	make_pair();
	margin_length = MARGIN;
}



