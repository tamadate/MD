//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------

void
MD::export_dump_recomb_CG(MD *md2) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int num1=int(vars->ions.size());
	int num=num1;
	if(flags->dump_gas==1) num+=int(vars->gases.size());
	if(flags->dump_2nd==1) {
		num+=int(md2->vars->ions.size());
		if(flags->dump_gas==1) num+=int(md2->vars->gases.size());
	}
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: ATOMS id type x y z\n", count, num);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {X=ion_r[0];Y=ion_r[1];Z=ion_r[2];}

	for (auto &a : vars->ions) fprintf(f, "%d %d %f %f %f\n", a.id, 100, ion_r[0]-X, ion_r[1]-Y, ion_r[2]-Z);
	if(flags->dump_gas) for (auto &a : vars->gases) fprintf(f, "%d %d %f %f %f\n", a.id, a.type, a.qx-X, a.qy-Y, a.qz-Z);
	if(flags->dump_2nd) {
		for (auto &a : md2->vars->ions) fprintf(f, "%d %d %f %f %f\n", a.id+num1, 101, md2->ion_r[0]-X, md2->ion_r[1]-Y, md2->ion_r[2]-Z);
		if(flags->dump_gas) for (auto &a : md2->vars->gases) fprintf(f, "%d %d %f %f %f\n", a.id+num1, a.type, a.qx-X, a.qy-Y, a.qz-Z);
	}
	fclose(f);
}



void
MD::export_dump(void) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int num1=int(vars->ions.size());
	int num=num1;
	if(flags->dump_gas==1) num+=int(vars->gases.size());

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z\n", count, num, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {X=ion_r[0];Y=ion_r[1];Z=ion_r[2];}

	for (auto &a : vars->ions) fprintf(f,"%d %d %f %f %f\n",a.id,a.type,a.qx-X,a.qy-Y,a.qz-Z);
	if(flags->dump_gas) for (auto &a : vars->gases) fprintf(f, "%d %d %f %f %f\n", a.id, a.type, a.qx-X, a.qy-Y, a.qz-Z);
	fclose(f);
}


void
MD::export_dump_recomb(MD *md2) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int num1=int(vars->ions.size());
	int num=num1;
	if(flags->dump_gas==1) num+=int(vars->gases.size());
	if(flags->dump_2nd==1) {
		num+=int(md2->vars->ions.size());
		if(flags->dump_gas==1) num+=int(md2->vars->gases.size());
	}
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: ATOMS id type x y z\n", count, num);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {X=ion_r[0];Y=ion_r[1];Z=ion_r[2];}

	for (auto &a : vars->ions) fprintf(f, "%d %d %f %f %f\n", a.id, a.type, a.qx-X, a.qy-Y, a.qz-Z);
	if(flags->dump_gas) for (auto &a : vars->gases) fprintf(f, "%d %d %f %f %f\n", a.id, a.type, a.qx-X, a.qy-Y, a.qz-Z);
	if(flags->dump_2nd) {
		for (auto &a : md2->vars->ions) fprintf(f, "%d %d %f %f %f\n", a.id+num1, a.type, a.qx-X, a.qy-Y, a.qz-Z);
		if(flags->dump_gas) for (auto &a : md2->vars->gases) fprintf(f, "%d %d %f %f %f\n", a.id+num1, a.type, a.qx-X, a.qy-Y, a.qz-Z);
	}
	fclose(f);
}

void
MD::export_dump_CG(void) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	int num1=int(vars->ions.size());
	int num=num1;
	if(flags->dump_gas==1) num+=int(vars->gases.size());

	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%e %e\n%e %e\n%e %e\nITEM: ATOMS id type x y z\n", count, num, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5, -d_size*0.5, d_size*0.5);

	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {X=ion_r[0];Y=ion_r[1];Z=ion_r[2];}

	for (auto &a : vars->ions) fprintf(f,"%d %d %f %f %f\n",a.id,a.type,a.qx,a.qy,a.qz);
	if(flags->dump_gas) for (auto &a : vars->gases) fprintf(f, "%d %d %f %f %f\n", a.id, a.type, a.qx-X, a.qy-Y, a.qz-Z);
	fclose(f);
}

void
MD::export_dump_rigid(MD *md2) {
	int count = vars->time;
	FILE*f=fopen(pp->dump_path, "a");
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n2\nITEM: ATOMS id type x y z\n", count);
	double X,Y,Z;
	X=Y=Z=0;
	if(flags->dump_fix==1) {X=ion_r[0];Y=ion_r[1];Z=ion_r[2];}

	fprintf(f, "%d %d %f %f %f\n", 1, 100, ion_r[0]-X, ion_r[1]-Y, ion_r[2]-Z);
	fprintf(f, "%d %d %f %f %f\n", 2, 101, md2->ion_r[0]-X, md2->ion_r[1]-Y, md2->ion_r[2]-Z);
	fclose(f);
}

