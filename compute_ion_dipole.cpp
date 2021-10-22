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
//	Ion *ions = vars->ions.data();
	const int gs = gas_in.size();
//	const int is = vars->ions.size();

	for (int k = 0; k < gs; k++) {
		int i = gas_in[k];
        double Ex,Ey,Ez,E0;
        Ex=Ey=Ez=0;
        for (auto &a : vars->ions) {
            double dx = gases[i].qx - a.qx;
            double dy = gases[i].qy - a.qy;
            double dz = gases[i].qz - a.qz;
            adjust_periodic(dx, dy, dz);
            Ex += a.charge/dx/dx;
            Ey += a.charge/dy/dy;
            Ez += a.charge/dz/dz;
        }
        E0=sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

    
        for (auto &a : vars->ions) {
            double dx = gases[i].qx - a.qx;
            double dy = gases[i].qy - a.qy;
            double dz = gases[i].qz - a.qz;
            double rsq = (dx * dx + dy * dy + dz * dz);
            double r2inv = 1/rsq;
            double rinv=sqrt(r2inv);
            double costh=(Ex*dx+Ey*dy+Ez*dz)*rinv/E0;
            double ion_dipole = -2.0*pp->alphagas*a.charge*r2inv*r2inv*qqrd2e*costh;
            gases[i].fx += ion_dipole*dx;
            gases[i].fy += ion_dipole*dy;
            gases[i].fz += ion_dipole*dz;
            a.fx -= ion_dipole*dx;
            a.fy -= ion_dipole*dy;
            a.fz -= ion_dipole*dz;
        }
	}
}

// charge devided by number of atoms in ions
void
MD::compute_ion_dipole2(void) {
	Gas *gases = vars->gases.data();
	const int gs = gas_in.size();
	const int is = vars->ions.size();
    
    for (int k = 0; k < gs; k++) {
        for (auto &a : vars->ions) {
            int i = gas_in[k];
            double dx = gases[i].qx - a.qx;
            double dy = gases[i].qy - a.qy;
            double dz = gases[i].qz - a.qz;
			adjust_periodic(dx, dy, dz);
			double rsq = (dx * dx + dy * dy + dz * dz);
            
            double r2inv = 1/rsq;
            double ion_dipole = -2*pp->alphagas*r2inv*r2inv*qqrd2e/is/is;
            double force_pair=ion_dipole*r2inv;
            gases[i].fx += force_pair * dx;
            gases[i].fy += force_pair * dy;
            gases[i].fz += force_pair * dz;
            a.fx -= force_pair * dx;
            a.fy -= force_pair * dy;
            a.fz -= force_pair * dz;
		}
	}
}


// Fitting of the kappa
void
MD::compute_ion_dipole3(void) {
	Gas *gases = vars->gases.data();
    Ion *ions = vars->ions.data();
	const int gs = gas_in.size();
    const int iion=30;

	for (int k = 0; k < gs; k++) {
        int i = gas_in[k];
        double dx = gases[i].qx - ions[iion].qx;
        double dy = gases[i].qy - ions[iion].qy;
        double dz = gases[i].qz - ions[iion].qz;
        adjust_periodic(dx, dy, dz);
        double rsq = (dx * dx + dy * dy + dz * dz);
        int type1=gases[i].type;
        int type2=ions[iion].type;
        if (rsq < CL2) {
            double r2inv = 1/rsq;
            double force_pair=-kappaIonDipole*r2inv*r2inv*r2inv*vars->pair_coeff[type1][type2][1]/double(pp->num_gas);
            gases[i].fx += force_pair * dx;
            gases[i].fy += force_pair * dy;
            gases[i].fz += force_pair * dz;
            ions[iion].fx -= force_pair * dx;
            ions[iion].fy -= force_pair * dy;
            ions[iion].fz -= force_pair * dz;
        }
	}
}


// Single charged site
void
MD::compute_ion_dipole4(void) {
	Gas *gases = vars->gases.data();
    Ion *ions = vars->ions.data();
	const int gs = gas_in.size();
    const int iion=30;
    
	for (int k = 0; k < gs; k++) {
        int i = gas_in[k];
        double dx = gases[i].qx - ions[iion].qx;
        double dy = gases[i].qy - ions[iion].qy;
        double dz = gases[i].qz - ions[iion].qz;
        adjust_periodic(dx, dy, dz);
        double rsq = (dx * dx + dy * dy + dz * dz);
        
        double r2inv = 1/rsq;
        double ion_dipole = -2*pp->alphagas*r2inv*r2inv*qqrd2e;
        double force_pair=ion_dipole*r2inv;

        gases[i].fx += force_pair * dx;
        gases[i].fy += force_pair * dy;
        gases[i].fz += force_pair * dz;
        ions[iion].fx -= force_pair * dx;
        ions[iion].fy -= force_pair * dy;
        ions[iion].fz -= force_pair * dz;
	}
}



