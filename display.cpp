#include "md.hpp"
/*########################################################################################

-----Display current properties-----

#######################################################################################*/


void
MD::display(int output_ONOFF){
	double gastemp = obs->gas_temperature(vars);
	double gastemp_total = obs->gas_total_temperature(vars, pp->num_gas); 
	double iontemp = obs->ion_temperature(vars);
	double virial=ters->compute_tersoff_virial(vars)/3.0/V*Cpress;
	double press=kb*vars->ions.size()*iontemp /(V*1e-30) + virial;
	double U = sw->energy ;
	double K = obs->ion_kinetic_energy(vars);
	std::cout << "----------------------TIME = " << vars->time/1000.0 << " ps-------------------------" << endl;
	printf("Kinetic energy = %1.2e		Temp(gas) = %f\nPotential energy = %1.2e		Temp(ion) = %f\nTotal energy = %1.2e			Temp(gas total) = %f\nPress(gas) = %f\n", K, gastemp, U, iontemp, K+U, gastemp_total, press/101300.0); 
	FILE*f=fopen("T-E.dat", "a");
	fprintf(f, "%e\t%e\n", iontemp, U);
	fclose(f);
}


