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


void
MD::computeN2(void) {
	Gas *gases = vars->gases.data();
    for(auto &i : gas_in){
        if(i<Nof_around_gas){
            int j=i+Nof_around_gas;
            double dx = gases[i].qx - gases[j].qx;
            double dy = gases[i].qy - gases[j].qy;
            double dz = gases[i].qz - gases[j].qz;
            double rsq = (dx * dx + dy * dy + dz * dz);
            double r = sqrt(rsq);
            double rk = 35000.0 * (r-1.098);
            double force_bond_harmonic = -2.0*rk/r;
            gases[i].fx += force_bond_harmonic * dx;
            gases[i].fy += force_bond_harmonic * dy;
            gases[i].fz += force_bond_harmonic * dz;
            gases[j].fx -= force_bond_harmonic * dx;
            gases[j].fy -= force_bond_harmonic * dy;
            gases[j].fz -= force_bond_harmonic * dz;
            
        
            /*double r2inv = 1/rsq;
            int type1=11;
            int type2=11;
			double r6inv = r2inv * r2inv * r2inv;
			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_pair = (force_lj)*r2inv;
			gases[i].fx += force_pair * dx;
			gases[i].fy += force_pair * dy;
			gases[i].fz += force_pair * dz;
			gases[j].fx -= force_pair * dx;
			gases[j].fy -= force_pair * dy;
			gases[j].fz -= force_pair * dz;*/
        }
	}
}


