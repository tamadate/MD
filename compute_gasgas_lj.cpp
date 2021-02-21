//------------------------------------------------------------------------
#include "md.hpp"
#include "variables.hpp"
#include "CalculationCond.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on gas-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_gasgas_lj(void) {
	Gas *gases = vars->gases.data();
	const int gs = gas_in.size();
	for(auto &a : pairs_gg){
		int i=a.i;
		int j=a.j;
		double dx = gases[i].qx - gases[j].qx;
		double dy = gases[i].qy - gases[j].qy;
		double dz = gases[i].qz - gases[j].qz;
		adjust_periodic(dx, dy, dz);
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=gases[i].type;
		int type2=gases[j].type;
		if (rsq < CL2) {
			double r6inv = r2inv * r2inv * r2inv;
			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_pair = (force_lj)*r2inv;
			gases[i].fx += force_pair * dx;
			gases[i].fy += force_pair * dy;
			gases[i].fz += force_pair * dz;
			gases[j].fx -= force_pair * dx;
			gases[j].fy -= force_pair * dy;
			gases[j].fz -= force_pair * dz;
		}	
	}
}


