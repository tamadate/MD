
//------------------------------------------------------------------------
#include "md.hpp"
#include "run.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*	
	- Precalculation for the recombination calculation			
	- MD simulation performed for (step_relax) steps.
	- Properties (pressure, temperature, energy etc...) are observed
	each (OBSERVE) steps.
*/
/////////////////////////////////////////////////////////////////////

void 
MD::run_recomb_pre(void) {
 	const int OBSERVE = 100000;
	for (int it=0; it < pp->step_relax; it++) {
		if (it%OBSERVE==0) {
			display(1);
			vars->time+=(OBSERVE*dt);
			export_dump();
		}
	   	verlet();
	}
}

/////////////////////////////////////////////////////////////////////
/*	
	- Re-precalculation for the next recombination simulation.
	- MD simulation performed for (step_repre) steps.
	- Properties (pressure, temperature, energy etc...) are observed
	each (OBSERVE) steps.
*/
/////////////////////////////////////////////////////////////////////

void 
MD::run_recomb_repre(void) {
 	const int OBSERVE = 100000;
	for (int it=0; it < pp->step_repre; it++) {
		if (it%OBSERVE==0) {
			display(1);
			vars->time+=(OBSERVE*dt);
		}
		verlet();
	}
}

/////////////////////////////////////////////////////////////////////
/*	
	- Recombination calculation
	- MD simulation is performed while collision or non-collision event
	is happen or after (runout) time.
	- Ions interaction between domain 1 and 2 is included, while 
	gas (domain1) - gas (domain2) interaction	is ignored. Then, 
	the density of the gas molecule is kept as unity value.
	- Properties (pressure, temperature, energy etc...) are observed
	each (OBSERVE) steps.
	- Md1 and Md2 maens each calculation domain.
	- return values i (-1, 0, 1) are (orbiting, non-collision, collision),
	respectively
*/
/////////////////////////////////////////////////////////////////////

int run_recomb(MD *md1, MD *md2, MD *MD1, MD *MD2) {
 	const int OBSERVE = 1000;
	int i=-1;
	int cg_f;
	double runout=100000000000;
	if (md1->flags->cg==0) cg_f=0;
	if (md1->flags->cg==1 || md1->flags->cg==2) cg_f=1;
	md1->rigid_dt=sqrt(md1->crsq)*10;
	while (i<0) {		// break loop when i=(0,1) <- (escape, collision)
		if (it%OBSERVE==0) {
			std::cout << "----------------------Domain 1-------------------------" << endl;
			md1->display(1);
			std::cout << "----------------------Domain 2-------------------------" << endl;
			md2->display(1);
			if(cg_f==0) md1->export_dump_recomb(md2);
			if(cg_f==1) {
				if (md1->flags->cg==1) md1->export_dump_recomb_CG(md2);
				if (md1->flags->cg==2) md1->export_dump_rigid(md2);
			}
			if(md1->flags->gyration==1 || md2->flags->gyration==1) md1->gyration_out(md2);
		}
   		if(cg_f==0)	{
			i=md1->verlet_recomb(md2);
			md1->vars->time+=dt;
			md2->vars->time+=dt;
		}
   		if(cg_f==1)	{
			if(md1->flags->cg==1) {
				i=md1->verlet_recomb_cg(md2);
				md1->vars->time+=dt;
				md2->vars->time+=dt;
			}
			if(md1->flags->cg==2) {
				i=md1->verlet_rigid(md2);
				md1->vars->time+=md1->rigid_dt;
				md2->vars->time+=md1->rigid_dt;
				md1->rigid_dt=sqrt(md1->crsq)*10;
			}
		}

   		if(cg_f==0 && md1->crsq>10000) {
			cg_f=1;
			md1->flex_to_rigid(md2, MD1, MD2);
		}
   		if(cg_f==1 && md1->crsq<10000) {
			cg_f=0; 
			MD1->run_recomb_repre();
			MD2->run_recomb_repre();
			md1->rigid_to_flex(md2, MD1, MD2);
		}

		if(md1->vars->time > runout) break;
	}
	return i;
}


/////////////////////////////////////////////////////////////////////
/*	
	- Calculation of ion properties after collision
	- MD calculation is performed for 10^5 steps and ion properties (velocity
	and structure) are logged.
	- ion velocities are saved on a directory (./velocity/)
	- ion strucutures are saved on a directory (./structure/)
*/
/////////////////////////////////////////////////////////////////////

void run_after(MD *md1, MD *md2, char* Ncalculation, int Ncollision) {
	sprintf(md1->filepath, "./velocity/%s_%d.dat", Ncalculation, Ncollision); 
	FILE*fv=fopen(md1->filepath, "w");
	sprintf(md1->filepath, "./structure/%s_%d.dump", Ncalculation, Ncollision); 
	FILE*fs=fopen(md1->filepath, "w");

	const int num1=int(md1->vars->ions.size());
	const int num2=int(md2->vars->ions.size());

	for (int it=0; it < 100000; it++) {
		if (it%100==0) {
			double v12=md1->ion_v[0]*md1->ion_v[0]+md1->ion_v[1]*md1->ion_v[1]+md1->ion_v[2]*md1->ion_v[2];
			double v22=md2->ion_v[0]*md2->ion_v[0]+md2->ion_v[1]*md2->ion_v[1]+md2->ion_v[2]*md2->ion_v[2];
			fprintf(fv, "%f\t%f\n", sqrt(v12), sqrt(v22));

			fprintf(fs, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: ATOMS id type x y z vx vy vz\n", int(md1->vars->time), num1+num2);
			for (auto &a : md1->vars->ions) fprintf(fs, "%d %d %f %f %f %f %f %f\n", a.id, a.type, a.qx-md1->ion_r[0], a.qy-md1->ion_r[1], a.qz-md1->ion_r[2], a.px, a.py, a.pz);
			for (auto &a : md2->vars->ions) fprintf(fs, "%d %d %f %f %f %f %f %f\n", a.id+num1, a.type, a.qx-md1->ion_r[0], a.qy-md1->ion_r[1], a.qz-md1->ion_r[2], a.px, a.py, a.pz);

		}
		int dammy=md1->verlet_recomb_after(md2);
	}
	fclose(fv);
	fclose(fs);
}










