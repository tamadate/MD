//------------------------------------------------------------------------
#include "md.hpp"
#include <random>
#include <algorithm>
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	

	- Initialization of recombinaiton simulation.
	- Here, the potision and velocites are delivered from Md to md
	- Two calculation domains are randomly rotated and set to initial 
	positions (0,0,0) and (0,delta,0)
	- Velocites are given as combinaiton value of thermal and electric 
	velocities (vth, vth+ve, vth)
	
*/
/////////////////////////////////////////////////////////////////////

int
MD::initial_recomb(MD *md2, MD *Md1, MD *Md2, Physical *pp12) {

//	take over positions/velocities of ion/gas (Md -> md)
	int gs=vars->gases.size();
	int is=vars->ions.size();
	for(int i=0;i<gs;i++){vars->gases[i]=Md1->vars->gases[i];}
	for(int i=0;i<is;i++){vars->ions[i]=Md1->vars->ions[i];}
	gs=md2->vars->gases.size();
	is=md2->vars->ions.size();
	for(int i=0;i<gs;i++){md2->vars->gases[i]=Md2->vars->gases[i];}
	for(int i=0;i<is;i++){md2->vars->ions[i]=Md2->vars->ions[i];}

//	set initial positions/velocities r=(0,0,0), vtrans=(0,0,0)
	for (auto &a : vars->ions) a.qx-=Md1->ion_r[0], a.qy-=Md1->ion_r[1], a.qz-=Md1->ion_r[2];
	for (auto &a : vars->gases) a.qx-=Md1->ion_r[0], a.qy-=Md1->ion_r[1], a.qz-=Md1->ion_r[2];
	for (auto &a : md2->vars->ions) a.qx-=Md2->ion_r[0], a.qy-=Md2->ion_r[1], a.qz-=Md2->ion_r[2];
	for (auto &a : md2->vars->gases) a.qx-=Md2->ion_r[0], a.qy-=Md2->ion_r[1], a.qz-=Md2->ion_r[2];
	analysis_ion();
	md2->analysis_ion();
	for (auto &a : vars->ions) a.px-=ion_v[0], a.py-=ion_v[1], a.pz-=ion_v[2];
	for (auto &a : md2->vars->ions) a.px-=md2->ion_v[0], a.py-=md2->ion_v[1], a.pz-=md2->ion_v[2];

//	Randomly rotating ions and change the center of domain 2 to (0,delta,0)
	random_device seed;
	double A,B,C,x,y,z;
	A=seed(),B=seed(),C=seed();
	for(auto &a : vars->ions) {x=a.qx,y=a.qy,z=a.qz; vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); x=a.px,y=a.py,z=a.pz; vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z);}
	for(auto &a : vars->gases) {x=a.qx,y=a.qy,z=a.qz; vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); x=a.px,y=a.py,z=a.pz; vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z);}
	A=seed(),B=seed(),C=seed();
	for(auto &a : md2->vars->ions) {x=a.qx,y=a.qy,z=a.qz; vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); x=a.px,y=a.py,z=a.pz; vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z); a.qy+=pp12->delta*1e10;}
	for(auto &a : md2->vars->gases) {x=a.qx,y=a.qy,z=a.qz; vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); x=a.px,y=a.py,z=a.pz; vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z); a.qy+=pp12->delta*1e10;}

//	add the thermal velocity as a transtrational velocity from MB distribution
//	add electrical velocity
	default_random_engine engine(seed());
	normal_distribution<> dist1(0.0, sqrt(kb*T/Md1->pp->m));	
	normal_distribution<> dist2(0.0, sqrt(kb*T/Md2->pp->m));
	double v01[3],v02[3];
	int flag=0;
	double ve1=Md1->pp->z*Md2->pp->z*e*e/4.0/M_PI/e0/pp12->delta/pp12->delta*Md1->pp->D/kb/T;
	double ve2=Md1->pp->z*Md2->pp->z*e*e/4.0/M_PI/e0/pp12->delta/pp12->delta*Md2->pp->D/kb/T;
	int n_kick=0;

	do{
		v01[1] = (dist1(engine)-ve1)*1e-5;
		v02[1] = (dist2(engine)+ve2)*1e-5;
		if(v02[1]-v01[1]<0.0){
			v01[0] = (dist1(engine))*1e-5;
			v01[2] = (dist1(engine))*1e-5;
			v02[0] = (dist2(engine))*1e-5;
			v02[2] = (dist2(engine))*1e-5;
			double vx=v01[0]-v02[0];
			double vy=v01[1]-v02[1];
			double vz=v01[2]-v02[2];
			pp12->v0=sqrt(vx*vx+vy*vy+vz*vz);
			pp12->th0=acos(vy/pp12->v0)/M_PI*180;
			double thc=pp12->thc(pp12 -> v0*1e5);
			n_kick++;
			if(thc*100000>pp12->th0) {
				flag++;
			}
		}
	} while (flag==0);

	for(auto &a : vars->ions) a.px+=v01[0], a.py+=v01[1], a.pz+=v01[2];
	for(auto &a : md2->vars->ions) a.px+=v02[0], a.py+=v02[1], a.pz+=v02[2];


	crsq=(pp12->delta*1e10)*(pp12->delta*1e10);

	double dy = pp12->delta*1e10;
	double force_pair = qqrd2e * pp->z * md2->pp->z / (dy*dy);
	ion_f[0] = 0;
	ion_f[1] = -force_pair;
	ion_f[2] = 0;
	md2->ion_f[0] = 0;
	md2->ion_f[1] = force_pair;
	md2->ion_f[2] = 0;


	for(auto &a : Md1->vflux) vflux.push_back(a);
	for(auto &a : Md2->vflux) md2->vflux.push_back(a);

	number=Md1->number;
	md2->number=Md2->number;
	analysis_ion();
	md2->analysis_ion();

	cout<<ion_v[0]<<" "<<ion_v[1]<<" "<<ion_v[2]<<" "<<endl;
	cout<<v01[0]<<" "<<v01[1]<<" "<<v01[2]<<" "<<endl;
	
	make_pair();
	if(flags->force_sw==1) sw->make_pair(vars);
	if(flags->force_ters==1) ters->make_pair(vars);
	md2->make_pair();
	if(md2->flags->force_sw==1) md2->sw->make_pair(vars);
	if(md2->flags->force_ters==1) md2->ters->make_pair(vars);

	margin_length = MARGIN;
	output_initial();
	del2=rmin2=md2->del2=md2->rmin2=(pp12->delta*1e10+10)*(pp12->delta*1e10+10);
	CD2=md2->CD2=pp12->CD*pp12->CD;


	setNVErecomb();
	md2->setNVErecomb();
	return n_kick;
}


