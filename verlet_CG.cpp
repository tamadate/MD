//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	- Compute a domain (NVE or NVT) for diffusion coefficient
	- v(t) -> v(t+dt/2)
	- Calculate velocity and potition of ion's center of mass.
	- Temperature control (velocity scaling, this case is NVT)
	- r(t) -> r(t+dt)
	- apply periodic boundary condition
	- Determine whethere update pair list or not
	- Calculate ion intraatmic interaction
	- Calculate ion-gas interatominc interaction 
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////
void
MD::verlet_cg(void) {
	velocity_calculation_cg();
	update_position_cg();	
	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	if(flags->nose_hoover_gas==1)	nosehoover_zeta_gas();
	check_pairlist();
	compute_cg();
	velocity_calculation_cg();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	if(flags->nose_hoover_gas==1)	nosehoover_gas();
}



/////////////////////////////////////////////////////////////////////
/*	
	- Compute two domains (semi-NVT) for recombination calculation
	- v(t) -> v(t+dt/2)
	- r(t) -> r(t+dt)
	- apply periodic boundary condition
	- Determine whethere update pair list or not
	- Calculate ion intraatomic interaction
	- Calculate ion-gas interatominc interaction 
	- Calculate domain-domain interaction and determine collision/non-
	collision
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////
int
MD::verlet_recomb_cg(MD *md2) {
	velocity_calculation_cg();
	md2->velocity_calculation_cg();

	update_position_cg();
	md2->update_position_cg();

	check_pairlist();
	md2->check_pairlist();

	compute_cg();
	md2->compute_cg();

	int judge=compute_domdom_rigid(md2);

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	md2->velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))

	return judge;
}




/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (half of a time step, dt/2)
*/
/////////////////////////////////////////////////////////////////////
void 
MD::velocity_calculation_cg(void) {
	Gas *gases = vars->gases.data();
	int gis=gas_in.size();
	double const Coeff=dt / 2.0 * 4.184 * 1e-4;
	ion_v[0] += ion_f[0] / pp->Mion *Coeff;
    ion_v[1] += ion_f[1] / pp->Mion *Coeff;
    ion_v[2] += ion_f[2] / pp->Mion *Coeff;
	for (int k=0;k<gis;k++){
		int i = gas_in[k];
		gases[i].px += gases[i].fx / gases[i].mass *Coeff;
		gases[i].py += gases[i].fy / gases[i].mass *Coeff;
		gases[i].pz += gases[i].fz / gases[i].mass *Coeff;
	}
}


/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (a time step, dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_position_cg(void) {
	Gas *gases = vars->gases.data();
	int gis=gas_in.size();

	ion_r[0] += ion_v[0]*dt;
    ion_r[1] += ion_v[1]*dt;
    ion_r[2] += ion_v[2]*dt;
	ion_f[0]=ion_f[1]=ion_f[2]=0.0;

	for (int k=0;k<gis;k++){
		for(int i = gas_in[k]; i < gas_in[k]+pp->num_gas; i++){
			gases[i].qx += gases[i].px * dt;
			gases[i].qy += gases[i].py * dt;
			gases[i].qz += gases[i].pz * dt;
			gases[i].fx=gases[i].fy=gases[i].fz=0.0;
		}
	}

}

/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-ion (Colulombic+sifted LJ)
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_cg(void) {
	Gas *gases = vars->gases.data();
	Ion *ions = vars->ions.data();
	const int gs = gas_in.size();
	for (int k = 0; k < gs; k++) {
		int i=gas_in[k];
		double dx=ion_r[0]-gases[i].qx;
		double dy=ion_r[1]-gases[i].qy;
		double dz=ion_r[2]-gases[i].qz;
		double rsq = (dx * dx + dy * dy + dz * dz);
		double r=sqrt(rsq);
		double shift_r=(r-pp->CG_r0);
		double r2inv=1/shift_r/shift_r;
		double r6inv = r2inv * r2inv * r2inv;
		int type1=gases[i].type;
		int type2=ions[0].type;
		double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		double force_pair = (force_lj)/r/shift_r;
		ion_f[0] += force_pair * dx;
		ion_f[1] += force_pair * dy;
		ion_f[2] += force_pair * dz;
		gases[i].fx -= force_pair * dx;
		gases[i].fy -= force_pair * dy;
		gases[i].fz -= force_pair * dz;
	}
}

