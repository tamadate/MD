//------------------------------------------------------------------------
#include "md.hpp"
#include "variables.hpp"
#include "CalculationCond.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_lj(void) {
	Gas *gases = vars->gases.data();
	Ion *ions = vars->ions.data();
	const int gs = gas_in.size();
	const int is = vars->ions.size();
	for(auto &a : pairs_gi){
		int i=a.i;
		int j=a.j;
		double dx = gases[i].qx - ions[j].qx;
		double dy = gases[i].qy - ions[j].qy;
		double dz = gases[i].qz - ions[j].qz;
		adjust_periodic(dx, dy, dz);
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=gases[i].type;
		int type2=ions[j].type;
		if (rsq < CL2) {
			double r6inv = r2inv * r2inv * r2inv;
			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_pair = (force_lj)*r2inv;
			gases[i].fx += force_pair * dx;
			gases[i].fy += force_pair * dy;
			gases[i].fz += force_pair * dz;
			ions[j].fx -= force_pair * dx;
			ions[j].fy -= force_pair * dy;
			ions[j].fz -= force_pair * dz;
		}	
	}

/*	for(auto &a : vars->ions){
		double force_pair = 5.55e-5*200.0;
		a.fx += force_pair;
	}*/
	
}


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-ion (Colulombic+LJ)
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_lj_coul(void) {
	Ion *ions = vars->ions.data();

	for (auto &a : vars->ion_pairs){
		int i=a.i;
		int j=a.j;
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r2inv = 1/rsq;
		int type1=ions[i].type;
		int type2=ions[j].type;
		double r6inv = r2inv * r2inv * r2inv;
		double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);	
		double force_coul = qqrd2e * ions[i].charge * ions[j].charge * sqrt(r2inv);
		double force_pair = (force_lj + force_coul)*r2inv;
		ions[i].fx += force_pair * dx;
		ions[i].fy += force_pair * dy;
		ions[i].fz += force_pair * dz;
		ions[j].fx -= force_pair * dx;
		ions[j].fy -= force_pair * dy;
		ions[j].fz -= force_pair * dz;
	}
}


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




