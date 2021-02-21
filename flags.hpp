#pragma once
#include "CalculationCond.hpp"

class FLAG{
public:
//NVT flags
	int velocity_scaling;
	int nose_hoover_ion;
	int nose_hoover_gas;
	int semi_NVT_gasgas;

//dump flags
	int dump_gas;
	int dump_fix;
	int dump_2nd;

//intra-atomic interaction in ion
	int force_ters;
	int force_sw;
	int force_coul;
	int intra_AMBER;

//inter-atomic interaction between ion-gas
	int force_lj;
	int force_ion_dipole;

//inter-atomic interaction between domain-domain
	int domdom_lj_coul;

//inter-atomic interaction between gas-gas
	int force_gasgas_lj;

//options
	int fix_cell_center;
	int gyration;
	int RDF;

//couase-grained
	int cg;


	FLAG(void);
	~FLAG(void);

private:

};

