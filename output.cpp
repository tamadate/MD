#include "md.hpp"
/*########################################################################################

-----Output-----

#######################################################################################*/

/**********************************initialization******************************************/
void
MD::output_initial(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number)); 
	f=fopen(filepath, "w");
	fclose(f);
}

void
MD::output(void){
	sprintf(filepath, "ion_%d_%d.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, ion_r[0], ion_r[1], ion_r[2], ion_v[0], ion_v[1], ion_v[2]);
	fclose(f);
}

void
MD::output_gas(void){
	sprintf(filepath, "gas_%d_%d.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%f %f %f %f %e %e %e\n", vars->time/1e6, gas_r[0], gas_r[1], gas_r[2], gas_v[0], gas_v[1], gas_v[2]);
	fclose(f);
}

void
MD::output_temp(double gastemp, double iontemp){
	sprintf(filepath, "%d_%d_relax.dat", int(T), int(calculation_number)); 
	FILE*f=fopen(filepath, "a");
	fprintf(f, "%f\t%f\t%f\n", vars->time/1e6, gastemp, iontemp);
	fclose(f);
}



