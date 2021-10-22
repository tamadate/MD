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
MD::v_verlet(void) {
	analysis_ion();
	if(flags->velocity_scaling==1)	velocity_scaling();
    
	update_position_v();	
	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	if(flags->nose_hoover_gas==1)	nosehoover_zeta_gas();
	check_pairlist();

	velocity_calculation_v();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	if(flags->nose_hoover_gas==1)	nosehoover_gas();
}


/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (half of a time step, dt/2)
*/
/////////////////////////////////////////////////////////////////////
void 
MD::velocity_calculation_v(void) {
	Gas *gases = vars->gases.data();
	double const Coeff=dt / 2.0 * 4.184 * 1e-4;
	for (auto &a : vars->ions) {
		a.px += a.fx / a.mass *Coeff;
	    a.py += a.fy / a.mass *Coeff;
	    a.pz += a.fz / a.mass *Coeff;
        a.fx=a.fy=a.fz=0;
	}
    for (auto &i : gas_in){
            gases[i].px += gases[i].fx / gases[i].mass *Coeff;
            gases[i].py += gases[i].fy / gases[i].mass *Coeff;
            gases[i].pz += gases[i].fz / gases[i].mass *Coeff;
            gases[i].fx=gases[i].fy=gases[i].fz=0;
    }
    compute_intra();
	compute_inter();
       
    for (auto &a : vars->ions) {
		a.px += a.fx / a.mass *Coeff;
	    a.py += a.fy / a.mass *Coeff;
	    a.pz += a.fz / a.mass *Coeff;
	}
    for (auto &i : gas_in){
        gases[i].px += gases[i].fx / gases[i].mass *Coeff;
        gases[i].py += gases[i].fy / gases[i].mass *Coeff;
        gases[i].pz += gases[i].fz / gases[i].mass *Coeff;
	}
    
    if(gastype==2)  update_velocity_constrained();
}

/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (a time step, dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_position_v(void) {
	Gas *gases = vars->gases.data();
    double const Coeff=dt * dt / 2.0 * 4.184 * 1e-4;
	for (auto &a : vars->ions) {
		a.qx += a.px * dt + a.fx / a.mass * Coeff;
		a.qy += a.py * dt + a.fy / a.mass * Coeff;
		a.qz += a.pz * dt + a.fz / a.mass * Coeff;
	}
    if(gastype==2)  update_position_constrained();
    else{
        for (auto &i : gas_in) {
            gases[i].qx += gases[i].px * dt + gases[i].fx / gases[i].mass * Coeff;
            gases[i].qy += gases[i].py * dt + gases[i].fy / gases[i].mass * Coeff;
            gases[i].qz += gases[i].pz * dt + gases[i].fz/ gases[i].mass * Coeff;
        }
	}
}

/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (a time step, dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_position_constrained(void) {
	Gas *gases = vars->gases.data();
    double const Coeff=dt * dt / 2.0 * 4.184 * 1e-4;
    for (auto &i : gas_in) {
        if(i<Nof_around_gas){
            int j=i+Nof_around_gas;
            double dx1=gases[j].qx-gases[i].qx;
            double dy1=gases[j].qy-gases[i].qy;
            double dz1=gases[j].qz-gases[i].qz;
            gases[i].qx += gases[i].px * dt + gases[i].fx / gases[i].mass * Coeff;
            gases[i].qy += gases[i].py * dt + gases[i].fy / gases[i].mass * Coeff;
            gases[i].qz += gases[i].pz * dt + gases[i].fz/ gases[i].mass * Coeff;
            gases[j].qx += gases[j].px * dt + gases[j].fx / gases[j].mass * Coeff;
            gases[j].qy += gases[j].py * dt + gases[j].fy / gases[j].mass * Coeff;
            gases[j].qz += gases[j].pz * dt + gases[j].fz / gases[j].mass * Coeff;
	
            double dx2=gases[j].qx-gases[i].qx;
            double dy2=gases[j].qy-gases[i].qy;
            double dz2=gases[j].qz-gases[i].qz;
            double r12=dx1*dx2+dy1*dy2+dz1*dz2;
            double r22=dx2*dx2+dy2*dy2+dz2*dz2;
            
            double ram=pp->Mgas*(r22-1.098*1.098)/(2*dt*dt*r12)/4.184e-4;
            gases[i].fx+=ram*dx2;
            gases[i].fy+=ram*dy2;
            gases[i].fz+=ram*dz2;
            gases[j].fx-=ram*dx2;
            gases[j].fy-=ram*dy2;
            gases[j].fz-=ram*dz2;
            
            gases[i].qx += ram * dx2 / gases[i].mass * Coeff;;
            gases[i].qy += ram * dy2 / gases[i].mass * Coeff;;
            gases[i].qz += ram * dz2 / gases[i].mass * Coeff;;
            gases[j].qx -= ram * dx2 / gases[j].mass * Coeff;;
            gases[j].qy -= ram * dy2 / gases[j].mass * Coeff;;
            gases[j].qz -= ram * dz2 / gases[j].mass * Coeff;;
        }
	}
}

void
MD::update_velocity_constrained(void) {
	Gas *gases = vars->gases.data();
    for (auto &i : gas_in) {
        if(i<Nof_around_gas){
            int j=i+Nof_around_gas;
            double const Coeff=dt / 2.0 * 4.184 * 1e-4;
            
            double dx1=gases[j].px-gases[i].px;
            double dy1=gases[j].py-gases[i].py;
            double dz1=gases[j].pz-gases[i].pz;
            double dx2=gases[j].qx-gases[i].qx;
            double dy2=gases[j].qy-gases[i].qy;
            double dz2=gases[j].qz-gases[i].qz;
            double r12=dx1*dx2+dy1*dy2+dz1*dz2;
            double r22=dx2*dx2+dy2*dy2+dz2*dz2;
            
            double ram=0.5*r12/(1.098*1.098);
            gases[i].px += ram * dx2;
            gases[i].py += ram * dy2;
            gases[i].pz += ram * dz2;
            gases[j].px -= ram * dx2;
            gases[j].py -= ram * dy2;
            gases[j].pz -= ram * dz2;
        }
	}
}



void
MD::v_verlet_eq(void) {
	analysis_ion();
	if(flags->velocity_scaling==1)	velocity_scaling();
	update_position_v();	
	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	if(flags->nose_hoover_gas==1)	nosehoover_zeta_gas();
	check_pairlist();

	velocity_calculation_veq();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	if(flags->nose_hoover_gas==1)	nosehoover_gas();
}

/////////////////////////////////////////////////////////////////////
/*	
	- Update velocity (half of a time step, dt/2)
*/
/////////////////////////////////////////////////////////////////////
void 
MD::velocity_calculation_veq(void) {
	Gas *gases = vars->gases.data();
	double const Coeff=dt / 2.0 * 4.184 * 1e-4;
	for (auto &a : vars->ions) {
		a.px += a.fx / a.mass *Coeff;
	    a.py += a.fy / a.mass *Coeff;
	    a.pz += a.fz / a.mass *Coeff;
        a.fx=a.fy=a.fz=0;
	}
    for (auto &i : gas_in){
            gases[i].px += gases[i].fx / gases[i].mass *Coeff;
            gases[i].py += gases[i].fy / gases[i].mass *Coeff;
            gases[i].pz += gases[i].fz / gases[i].mass *Coeff;
            gases[i].fx=gases[i].fy=gases[i].fz=0;
    }
    
    compute_intra();
       
    for (auto &a : vars->ions) {
		a.px += a.fx / a.mass *Coeff;
	    a.py += a.fy / a.mass *Coeff;
	    a.pz += a.fz / a.mass *Coeff;
	}
    for (auto &i : gas_in){
        gases[i].px += gases[i].fx / gases[i].mass *Coeff;
        gases[i].py += gases[i].fy / gases[i].mass *Coeff;
        gases[i].pz += gases[i].fz / gases[i].mass *Coeff;
	}
    
    if(gastype==2)  update_velocity_constrained();
}
