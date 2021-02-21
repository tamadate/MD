//------------------------------------------------------------------------
#include "md.hpp"

/**********************************constructor******************************************/
FLAG::FLAG(void) {
	velocity_scaling=0;
	nose_hoover_ion=0;
	nose_hoover_gas=0;
	semi_NVT_gasgas=0;
	dump_gas=1;
	dump_fix=1;
	fix_cell_center=0;
	dump_2nd=1;
	gyration=0;
	RDF=0;

//	ion intratomic interaction
	force_ters=0;
	force_sw=0;
	force_coul=0;
	intra_AMBER==0;

//	ion-gas interaction
	force_lj=0;
	force_ion_dipole=0;

//	domain-domain interaction
	domdom_lj_coul=0;

//	gas-gas interaction
	force_gasgas_lj=0;

//couase-grain
	cg=0;


}

/**********************************destractor******************************************/
FLAG::~FLAG(void) {

}


