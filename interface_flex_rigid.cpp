//------------------------------------------------------------------------
#include "md.hpp"
#include <random>
#include <algorithm>
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	Interface from rigid(CG) MD to all-atom MD
*/
/////////////////////////////////////////////////////////////////////
void
MD::rigid_to_flex(MD *md2, MD *Md1, MD *Md2) {

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
	Md1->analysis_ion();
	Md2->analysis_ion();
	for (auto &a : vars->ions) a.px-=Md1->ion_v[0], a.py-=Md1->ion_v[1], a.pz-=Md1->ion_v[2];
	for (auto &a : md2->vars->ions) a.px-=Md2->ion_v[0], a.py-=Md2->ion_v[1], a.pz-=Md2->ion_v[2];


//	Randomly rotating ions and change the center of domain 2 to (0,delta,0)
	random_device seed;
	double A,B,C,x,y,z;
	A=seed(),B=seed(),C=seed();
	for(auto &a : vars->ions) {
		x=a.qx,y=a.qy,z=a.qz; 
		vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); 
		x=a.px,y=a.py,z=a.pz; 
		vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z);
		a.qx+=ion_r[0];
		a.qy+=ion_r[1];
		a.qz+=ion_r[2];
	}
	for(auto &a : vars->gases) {
		x=a.qx,y=a.qy,z=a.qz; 
		vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); 
		x=a.px,y=a.py,z=a.pz; 
		vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z);
		a.qx+=ion_r[0];
		a.qy+=ion_r[1];
		a.qz+=ion_r[2];
	}
	A=seed(),B=seed(),C=seed();
	for(auto &a : md2->vars->ions) {
		x=a.qx,y=a.qy,z=a.qz; 
		vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); 
		x=a.px,y=a.py,z=a.pz; 
		vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z); 
		a.qx+=md2->ion_r[0];
		a.qy+=md2->ion_r[1];
		a.qz+=md2->ion_r[2];
	}
	for(auto &a : md2->vars->gases) {
		x=a.qx,y=a.qy,z=a.qz; 
		vars->ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z); 
		x=a.px,y=a.py,z=a.pz; 
		vars->ROTATION(a.px,a.py,a.pz,A,B,C,x,y,z); 
		a.qx+=md2->ion_r[0];
		a.qy+=md2->ion_r[1];
		a.qz+=md2->ion_r[2];
	}

//	add the thermal velocity as a transtrational velocity from MB distribution
//	add electrical velocity

	for(auto &a : vars->ions) a.px+=ion_v[0], a.py+=ion_v[1], a.pz+=ion_v[2];
	for(auto &a : md2->vars->ions) a.px+=md2->ion_v[0], a.py+=md2->ion_v[1], a.pz+=md2->ion_v[2];


	make_pair();
	if(flags->force_sw==1) sw->make_pair(vars);
	if(flags->force_ters==1) ters->make_pair(vars);
	md2->make_pair();
	if(md2->flags->force_sw==1) md2->sw->make_pair(vars);
	if(md2->flags->force_ters==1) md2->ters->make_pair(vars);

	margin_length = MARGIN;

	compute_intra();
	md2->compute_intra();

	compute_inter();
	md2->compute_inter();
	int judge=compute_domdom(md2);

}



/////////////////////////////////////////////////////////////////////
/*	
	Interface from all-atom MD to rigid(CG) MD
*/
/////////////////////////////////////////////////////////////////////
void
MD::flex_to_rigid(MD *md2, MD *Md1, MD *Md2) {
	ion_f[0]=ion_f[1]=ion_f[2]=0.0;
	md2->ion_f[0]=md2->ion_f[1]=md2->ion_f[2]=0.0;
	int dammy=compute_domdom_rigid(md2);
}


