//------------------------------------------------------------------------
#include "md.hpp"
#include "variables.hpp"
#include "CalculationCond.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the interaction acting on ion-ion (LJ+Coulombic)
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_domdom_lj_coul(MD *md2) {
	for (auto &a : vars->ions) {
		for (auto &a2 : md2->vars->ions) {
			double dx = a2.qx - a.qx;
			double dy = a2.qy - a.qy;
			double dz = a2.qz - a.qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			int type1=a2.type;
			int type2=a.type;
			double r6inv = r2inv * r2inv * r2inv;

			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_coul = qqrd2e * a.charge * a2.charge * sqrt(r2inv);
			double force_pair = (force_lj + force_coul)*r2inv;
			a2.fx += force_pair * dx;
			a2.fy += force_pair * dy;
			a2.fz += force_pair * dz;
			a.fx -= force_pair * dx;
			a.fy -= force_pair * dy;
			a.fz -= force_pair * dz;
		}
	}
}


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the interaction acting on ion-ion (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_domdom_lj(MD *md2) {
	for (auto &a : vars->ions) {
		for (auto &a2 : md2->vars->ions) {
			double dx = a2.qx - a.qx;
			double dy = a2.qy - a.qy;
			double dz = a2.qz - a.qz;
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			int type1=a2.type;
			int type2=a.type;
			double r6inv = r2inv * r2inv * r2inv;

			double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
			double force_pair = (force_lj)*r2inv;
			a2.fx += force_pair * dx;
			a2.fy += force_pair * dy;
			a2.fz += force_pair * dz;
			a.fx -= force_pair * dx;
			a.fy -= force_pair * dy;
			a.fz -= force_pair * dz;
		}
	}
}


