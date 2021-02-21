//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------



/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the ion induced dipole force working between ion and gas 
	molecules. Because ion have a charge, ion induce the gas dipole.
*/
/////////////////////////////////////////////////////////////////////

void
MD::compute_ion_dipole(void) {
	Gas *gases = vars->gases.data();
	Ion *ions = vars->ions.data();
	const int gs = gas_in.size();
	const int is = vars->ions.size();

	for (int k = 0; k < gs; k++) {
		for(int i = gas_in[k]; i < gas_in[k]+pp->num_gas; i++){
			double dx = gases[i].qx - ion_r[0];
			double dy = gases[i].qy - ion_r[1];
			double dz = gases[i].qz - ion_r[2];
			adjust_periodic(dx, dy, dz);
			double rsq = (dx * dx + dy * dy + dz * dz);
			if (rsq < CL2) {
				double r2inv = 1/rsq;
				double ion_dipole = -4.0*pp->alphagas/2.0*r2inv*r2inv*qqrd2e;
				double force_pair=ion_dipole*r2inv;
				gases[i].fx += force_pair * dx;
				gases[i].fy += force_pair * dy;
				gases[i].fz += force_pair * dz;
				force_pair/=double(is);
				for (auto &a : vars->ions) {
					a.fx -= force_pair * dx;
					a.fy -= force_pair * dy;
					a.fz -= force_pair * dz;
				}
			}
		}
	}
}


