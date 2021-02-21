//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	- Calculate the intra-atomic interaction.
*/
/////////////////////////////////////////////////////////////////////
void
MD::compute_intra(void) {
	if(flags->force_sw==1) sw->compute_sw(vars);
	if(flags->force_ters==1) ters->compute_tersoff(vars);
	if(flags->force_coul==1) compute_coul();
	if(flags->intra_AMBER==1){
		compute_lj_coul();
		bond_harmonic();
		angle_harmonic();
		dihedral_fourier();
	}
}
