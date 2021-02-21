#include "md.hpp"
/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/

/**********************************Force calculation******************************************/
void
MD::bond_harmonic(void) {
	Ion *ions = vars->ions.data();
/*intra-molecular interaction (bond)*/
	Bond_type *btypes = vars->btypes.data();
	for (auto &b : vars-> bonds) {
		int i=b.atom1, j=b.atom2, type=(b.type);
		double dx = ions[i].qx - ions[j].qx;
		double dy = ions[i].qy - ions[j].qy;
		double dz = ions[i].qz - ions[j].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r = sqrt(rsq);
		double rk = btypes[type].coeff1 * (r-btypes[type].coeff2);
		double force_bond_harmonic;
		force_bond_harmonic = -2.0*rk/r;
		ions[i].fx += force_bond_harmonic * dx;
		ions[i].fy += force_bond_harmonic * dy;
		ions[i].fz += force_bond_harmonic * dz;
		ions[j].fx -= force_bond_harmonic * dx;
		ions[j].fy -= force_bond_harmonic * dy;
		ions[j].fz -= force_bond_harmonic * dz;
	}
}
