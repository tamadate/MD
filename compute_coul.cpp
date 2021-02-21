#include "md.hpp"
#include "variables.hpp"
#include "CalculationCond.hpp"
/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the ion induced dipole force working between ion and gas 
	molecules. Because ion have a charge, ion induce the gas dipole.
*/
/////////////////////////////////////////////////////////////////////

void
MD::compute_coul(void) {
	Ion *ions = vars->ions.data();
	const int is = vars->ions.size();

	for (int i = 0; i < is-1; i++) {
		for(int j = i+1; j < is; j++){
			double dx = ions[i].qx - ions[j].qx;
			double dy = ions[i].qy - ions[j].qy;
			double dz = ions[i].qz - ions[j].qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double force_coul = qqrd2e * ions[j].charge * ions[i].charge * sqrt(r2inv);
			double force_pair = (force_coul)*r2inv;
			ions[i].fx += force_pair * dx;
			ions[i].fy += force_pair * dy;
			ions[i].fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
		}
	}
}


