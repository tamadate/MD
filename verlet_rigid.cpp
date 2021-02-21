//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	**Caution** this simulation can apply only on without gas calculation

	- Compute two domains (semi-NVT)
	- v(t) -> v(t+dt/2)
	- r(t) -> r(t+dt)
	- Calculate domain-domain interaction and determine collision/non-
	collision
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////
int
MD::verlet_rigid(MD *md2) {
	double const Coeff=rigid_dt / 2.0 * 4.184 * 1e-4;
	ion_v[0] += ion_f[0] / pp->Mion *Coeff;
    ion_v[1] += ion_f[1] / pp->Mion *Coeff;
    ion_v[2] += ion_f[2] / pp->Mion *Coeff;
	md2->ion_v[0] += md2->ion_f[0] / md2->pp->Mion *Coeff;
    md2->ion_v[1] += md2->ion_f[1] / md2->pp->Mion *Coeff;
    md2->ion_v[2] += md2->ion_f[2] / md2->pp->Mion *Coeff;

	ion_r[0] += ion_v[0]*rigid_dt;
    ion_r[1] += ion_v[1]*rigid_dt;
    ion_r[2] += ion_v[2]*rigid_dt;
	md2->ion_r[0] += md2->ion_v[0]*rigid_dt;
    md2->ion_r[1] += md2->ion_v[1]*rigid_dt;
    md2->ion_r[2] += md2->ion_v[2]*rigid_dt;
	ion_f[0]=ion_f[1]=ion_f[2]=0.0;
	md2->ion_f[0]=md2->ion_f[1]=md2->ion_f[2]=0.0;

	int judge=compute_domdom_rigid(md2);

	ion_v[0] += ion_f[0] / pp->Mion *Coeff;
    ion_v[1] += ion_f[1] / pp->Mion *Coeff;
    ion_v[2] += ion_f[2] / pp->Mion *Coeff;
	md2->ion_v[0] += md2->ion_f[0] / md2->pp->Mion *Coeff;
    md2->ion_v[1] += md2->ion_f[1] / md2->pp->Mion *Coeff;
    md2->ion_v[2] += md2->ion_f[2] / md2->pp->Mion *Coeff;	//	v(t+rigid_dt/2) -> v(t+rigid_dt) using F(x(t+rigid_dt))

	return judge;
}



/////////////////////////////////////////////////////////////////////
/*	
	- Calculate force working on ion-ion (Colulombic)
*/
/////////////////////////////////////////////////////////////////////
int
MD::compute_domdom_rigid(MD *md2) {
	double dx=ion_r[0]-md2->ion_r[0];
	double dy=ion_r[1]-md2->ion_r[1];
	double dz=ion_r[2]-md2->ion_r[2];
	double rsq = (dx * dx + dy * dy + dz * dz);
	double r2inv = 1/rsq;
	double force_pair = qqrd2e * pp->z * md2->pp->z * sqrt(r2inv) * r2inv;
	ion_f[0] += force_pair * dx;
	ion_f[1] += force_pair * dy;
	ion_f[2] += force_pair * dz;
	md2->ion_f[0] -= force_pair * dx;
	md2->ion_f[1] -= force_pair * dy;
	md2->ion_f[2] -= force_pair * dz;

	int judge=-1;
	crsq=rsq;
	if (crsq < rmin2) rmin2=crsq;
	if (crsq > del2) judge=0;
	if (crsq < CD2) judge=1;
	return judge;
}


