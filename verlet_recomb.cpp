//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*	
	- Compute two domains (semi-NVT)
	- v(t) -> v(t+dt/2)
	- r(t) -> r(t+dt)
	- apply periodic boundary condition
	- Determine whethere update pair list or not
	- Calculate ion intraatomic interaction
	- Calculate ion-gas interatominc interaction 
	- Calculate domain-domain interaction and determine collision/non-
	collision
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////
int
MD::verlet_recomb(MD *md2) {
	velocity_calculation();
	md2->velocity_calculation();

	analysis_ion();
	md2->analysis_ion();
	update_position();
	md2->update_position();
	periodic();
	md2->periodic();

	check_pairlist();
	md2->check_pairlist();

	compute_intra();
	md2->compute_intra();

	compute_inter();
	md2->compute_inter();
	int judge=compute_domdom(md2);

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	md2->velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))

	return judge;
}


/////////////////////////////////////////////////////////////////////
/*	
	- Compute two domains (semi-NVT)
	- v(t) -> v(t+dt/2)
	- r(t) -> r(t+dt)
	- apply periodic boundary condition
	- Determine whethere update pair list or not
	- Calculate ion intraatomic interaction
	- Calculate ion-gas interatominc interaction 
	- Calculate domain-domain interaction and determine collision/non-
	collision
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////

int
MD::verlet_recomb_after(MD *md2) {
	velocity_calculation();
	md2->velocity_calculation();

	analysis_ion();
	md2->analysis_ion();
	update_position();
	md2->update_position();
	periodic();
	md2->periodic();

	check_pairlist();
	md2->check_pairlist();

	compute_intra();
	md2->compute_intra();

	compute_inter();
	md2->compute_inter();
	int judge=compute_domdom_after(md2);

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))
	md2->velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))

	return judge;
}

