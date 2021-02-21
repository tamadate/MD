//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	- Compute a domain (NVE or NVT)
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
MD::verlet(void) {
	velocity_calculation();
	analysis_ion();
	if(flags->fix_cell_center==1) fix_cell_center(); // off
	if(flags->velocity_scaling==1)	velocity_scaling();

	update_position();	
	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	if(flags->nose_hoover_gas==1)	nosehoover_zeta_gas();
	if(flags->semi_NVT_gasgas==0)	periodic();
	if(flags->semi_NVT_gasgas==1)	periodic_semi();
	check_pairlist();
	compute_intra();
	compute_inter();

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	if(flags->nose_hoover_gas==1)	nosehoover_gas();
}


/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (half of a time step, dt/2)
*/
/////////////////////////////////////////////////////////////////////
void 
MD::velocity_calculation(void) {
	Gas *gases = vars->gases.data();
	int gis=gas_in.size();
	double const Coeff=dt / 2.0 * 4.184 * 1e-4;
	for (auto &a : vars->ions) {
		a.px += a.fx / a.mass *Coeff;
	    a.py += a.fy / a.mass *Coeff;
	    a.pz += a.fz / a.mass *Coeff;
	}
	for (int k=0;k<gis;k++){
		for(int i = gas_in[k]; i < gas_in[k]+pp->num_gas; i++){
			gases[i].px += gases[i].fx / gases[i].mass *Coeff;
			gases[i].py += gases[i].fy / gases[i].mass *Coeff;
			gases[i].pz += gases[i].fz / gases[i].mass *Coeff;
		}
	}
}

/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (a time step, dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_position(void) {
	Gas *gases = vars->gases.data();
	int gis=gas_in.size();
	for (auto &a : vars->ions) {
		a.qx += a.px * dt;
		a.qy += a.py * dt;
		a.qz += a.pz * dt;
		a.fx=a.fy=a.fz=0.0;
	}
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
	- Fix center of domain (v-v_center)
*/
/////////////////////////////////////////////////////////////////////
void
MD::fix_cell_center(void) {
	for (auto &a : vars->ions) {
		a.px -= ion_v[0];
		a.py -= ion_v[1];
		a.pz -= ion_v[2];
	}
}

