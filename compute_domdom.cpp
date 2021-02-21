//------------------------------------------------------------------------
#include "md.hpp"
#include "variables.hpp"
#include "CalculationCond.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the interaction acting on inter-domain. Normally L-J 
	and	Coulombic potentials based ion-ion interactions are calculated.
*/
/////////////////////////////////////////////////////////////////////
int
MD::compute_domdom(MD *md2) {
	int judge=-1;
	compute_domdom_lj_coul(md2);
	analysis_ion();
	double dx=ion_r[0]-md2->ion_r[0];
	double dy=ion_r[1]-md2->ion_r[1];
	double dz=ion_r[2]-md2->ion_r[2];
	crsq=(dx * dx + dy * dy + dz * dz);
	if (crsq < rmin2) rmin2=crsq;
	if (crsq > del2) judge=0;
	if (crsq < CD2) judge=1;
	return judge;
}


int
MD::compute_domdom_after(MD *md2) {
	int judge=-1;
	compute_domdom_lj(md2);
	analysis_ion();
	double dx=ion_r[0]-md2->ion_r[0];
	double dy=ion_r[1]-md2->ion_r[1];
	double dz=ion_r[2]-md2->ion_r[2];
	crsq=(dx * dx + dy * dy + dz * dz);
	if (crsq < rmin2) rmin2=crsq;
	if (crsq > del2) judge=0;
	if (crsq < CD2) judge=1;
	return judge;
}



