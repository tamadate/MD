//------------------------------------------------------------------------
#include "md.hpp"
#include "variables.hpp"
#include "CalculationCond.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the interatmic interaction working on ion-gas.
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_inter(void) {
	if(flags->force_lj==1)	compute_lj();
	if(flags->force_ion_dipole==1)	compute_ion_dipole();
	if(flags->force_gasgas_lj==1)	compute_gasgas_lj();
}


