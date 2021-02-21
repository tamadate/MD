//------------------------------------------------------------------------
#include "md.hpp"
#include "variables.hpp"
#include "CalculationCond.hpp"
//------------------------------------------------------------------------



/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-ion (Colulombic+LJ)
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_lj_coul(void) {
	Gas *gases = vars->gases.data();
	Ion *ions = vars->ions.data();
	const int gs = gas_in.size();
	const int is = vars->ions.size();

	for (auto &a : pairs){
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


