//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	constructor
*/
/////////////////////////////////////////////////////////////////////
MD::MD(char* condfile, char* file) {
	vars = new Variables();
	obs = new Observer();
	pp = new Physical();
	ters = new Tersoff();
	sw = new SW();
	flags = new FLAG();
	pp->PhysicalProp_set(condfile, file, flags);
	vars->read_initial(file);
	vars->set_initial_velocity(pp);
	initialization_gas(vars);
	analysis_ion();
	make_pair();
	margin_length = MARGIN;


}

/////////////////////////////////////////////////////////////////////
/*	
	destructor
*/
/////////////////////////////////////////////////////////////////////
MD::~MD(void) {
	delete vars;
	delete obs;
	delete pp;
	delete ters;
	delete flags;
}


