#include "md.hpp"
#include <random>
#include <algorithm>
/*########################################################################################

-----Output-----

#######################################################################################*/

/**********************************initialization******************************************/
void
MD::velocity_distribution(void){
	double vmax=5.0*pp->cgas*1e-5;
	double dv=vmax/MM;	
	int j;
	for(auto &a : vars->gases){
		double v2=a.px*a.px+a.py*a.py+a.pz*a.pz;
		double vtotal=sqrt(v2);
		j=vtotal/dv;
		disv[j]++;
		if(a.px>0){j=a.px/dv;disvx_pos[j]++;}
		else{j=-a.px/dv;disvx_neg[j]++;}
		if(a.py>0){j=a.py/dv;disvy_pos[j]++;}
		else{j=-a.py/dv;disvy_neg[j]++;}
		if(a.pz>0){j=a.pz/dv;disvz_pos[j]++;}
		else{j=-a.pz/dv;disvz_neg[j]++;}
	}
}

void
MD::update_distribution(void){
	double vmax=5.0*pp->cgas;
	double dv=vmax/MM;	
	double TIME=vars->time*1e-6;
	sprintf(filepath, "gas_dist_%d.dat", int(TIME)); 
	FILE*f=fopen(filepath, "w");
	for(int i=0; i<MM; i++) {
		fprintf(f, "%f\t%f\n", dv*(i+0.5), disv[i]/Nof_around_gas/2000/dv);
		disv[i]=0;
	}
	fclose(f);
	sprintf(filepath, "gas_dist_xyz_%d.dat", int(TIME)); 
	f=fopen(filepath, "w");
	for(int i=0; i<MM; i++) {
		fprintf(f, "%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n",  dv*(i+0.5), disvx_pos[i]/Nof_around_gas/dv/UPDATE_DISTRIBUTION*VELOCITY_DISTRIBUTION, disvy_pos[i]/Nof_around_gas/dv/UPDATE_DISTRIBUTION*VELOCITY_DISTRIBUTION, disvz_pos[i]/Nof_around_gas/dv/UPDATE_DISTRIBUTION*VELOCITY_DISTRIBUTION, -dv*(i+0.5), disvx_neg[i]/Nof_around_gas/dv/UPDATE_DISTRIBUTION*VELOCITY_DISTRIBUTION, disvy_neg[i]/Nof_around_gas/dv/UPDATE_DISTRIBUTION*VELOCITY_DISTRIBUTION, disvz_neg[i]/Nof_around_gas/dv/UPDATE_DISTRIBUTION*VELOCITY_DISTRIBUTION);
		disvx_pos[i]=disvy_pos[i]=disvz_pos[i]=disvx_neg[i]=disvy_neg[i]=disvz_neg[i]=0;
	}
	fclose(f);
}
