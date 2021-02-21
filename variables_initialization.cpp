//------------------------------------------------------------------------
#include "variables.hpp"
//------------------------------------------------------------------------
void 
Variables::read_initial(char* infile) {
	ifstream stream(infile);
	string str;
	int iflag=0;
	int num_atoms=0;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type name mass coeff1 coeff2") {iflag=1; continue;}
		if (str=="bond type name coeff1 coeff2") {iflag=2; continue;}
		if (str=="angle type name coeff1 coeff2") {iflag=3; continue;}
		if (str=="dihedral type name coeff1 coeff2 coeff3 coeff4") {iflag=4; continue;}
		if (str=="atoms") {iflag=5; continue;}
		if (str=="bonds") {iflag=6; continue;}
		if (str=="angles") {iflag=7; continue;}
		if (str=="dihedrals") {iflag=8; continue;}
		if (str=="molecules") {iflag=9; continue;}
	    string tmp;
    	istringstream stream(str);		
		if (iflag==1) {
			int loop=0;
			int type;
			double mass, coeff1, coeff2;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) type=stoi(tmp);
				if (loop==2) mass=stod(tmp);
				if (loop==3) coeff1=stod(tmp);
				if (loop==4) coeff2=stod(tmp);
				loop++;
			}
			add_atype(type, mass, coeff1, coeff2);
			num_atoms++;
		}
		if (iflag==2) {
			int loop=0;
			int type;
			double coeff1, coeff2;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) type=stoi(tmp);
				if (loop==2) coeff1=stod(tmp);
				if (loop==3) coeff2=stod(tmp);
				loop++;
			}
			add_btype(type, coeff1, coeff2);
		}
		if (iflag==3) {
			int loop=0;
			int type;
			double coeff1, coeff2;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) type=stoi(tmp);
				if (loop==2) coeff1=stod(tmp);
				if (loop==3) coeff2=stod(tmp)/180.0*M_PI;
				loop++;
			}
			add_ctype(type, coeff1, coeff2);
		}
		if (iflag==4) {
			int loop=0;
		 	int type, multi;
			double coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8, coeff9, coeff10,  coeff11, coeff12, coeff13, coeff14, coeff15;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) type=stoi(tmp);
				if (loop==2) multi=stod(tmp);
				if (loop==2) coeff1=stod(tmp);
				if (loop==3) coeff2=stod(tmp);
				if (loop==4) {
					coeff3=stod(tmp);
					coeff4=cos(stod(tmp)/180.0*M_PI);
					coeff5=sin(stod(tmp)/180.0*M_PI);
				}
				if (loop==5) coeff6=stod(tmp);
				if (loop==6) coeff7=stod(tmp);
				if (loop==7) {
					coeff8=stod(tmp);
					coeff9=cos(stod(tmp)/180.0*M_PI);
					coeff10=sin(stod(tmp)/180.0*M_PI);
				}
				if (loop==8) coeff11=stod(tmp);
				if (loop==9) coeff12=stod(tmp);
				if (loop==10) {
					coeff13=stod(tmp);
					coeff14=cos(stod(tmp)/180.0*M_PI);
					coeff15=sin(stod(tmp)/180.0*M_PI);
				}
				loop++;
			}
			add_dtype(type, multi, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8, coeff9, coeff10, coeff11, coeff12, coeff13, coeff14, coeff15);
		}
		if (iflag==5) {
			int loop=0;
			int type, id;
			double charge, x, y, z, vx, vy, vz, mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) id=stoi(tmp);
				if (loop==1) type=stoi(tmp);
				if (loop==2) charge=stod(tmp);
				if (loop==3) x=stod(tmp);
				if (loop==4) y=stod(tmp);
				if (loop==5) z=stod(tmp);
				if (loop==6) vx=stod(tmp);
				if (loop==7) vy=stod(tmp);
				if (loop==8) vz=stod(tmp);
				loop++;
			}
			for (auto &at : atypes) {
				if(at.type==type) mass=at.mass;
			}
			add_ions(id,type,x,y,z,vx,vy,vz,0,0,0,charge,mass);
		}
		if (iflag==6) {
			int loop=0;
			int atom1, atom2, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) atom1=stoi(tmp);
				if (loop==1) atom2=stoi(tmp);
				if (loop==2) type=stoi(tmp);
				loop++;
			}
			add_bonds(atom1-1,atom2-1,type-1);
		}
		if (iflag==7) {
			int loop=0;
			int atom1, atom2, atom3, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) atom1=stoi(tmp);
				if (loop==1) atom2=stoi(tmp);
				if (loop==2) atom3=stoi(tmp);
				if (loop==3) type=stoi(tmp);
				loop++;
			}
			add_angles(atom1-1,atom2-1,atom3-1,type-1);
		}
		if (iflag==8) {
			int loop=0;
			int atom1, atom2, atom3, atom4, type;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) atom1=stoi(tmp);
				if (loop==1) atom2=stoi(tmp);
				if (loop==2) atom3=stoi(tmp);
				if (loop==3) atom4=stoi(tmp);
				if (loop==4) type=stoi(tmp);
				loop++;
			}
			add_dihedrals(atom1-1,atom2-1,atom3-1,atom4-1,type-1);
		}
		if (iflag==9) {
			int atom;
			atom=stoi(str);
			molecules.push_back(atom);
		}
	}
	pair_coeff.resize(num_atoms+1);
	for (int i=0;i<num_atoms+1;i++){
		pair_coeff[i].resize(num_atoms+1);
		for (int j=0;j<num_atoms+1;j++){
			pair_coeff[i][j].resize(2);
		}
	}
	for (int i=0;i<num_atoms;i++){
		for (int j=0;j<num_atoms;j++){
			double epu, sigma;
			epu=sqrt(atypes[i].coeff1*atypes[j].coeff1);
			sigma=(atypes[i].coeff2+atypes[j].coeff2)*0.5;
			pair_coeff[atypes[i].type][atypes[j].type][0]=48 * epu*pow(sigma,12.0);
			pair_coeff[atypes[i].type][atypes[j].type][1]=24 * epu*pow(sigma,6.0);
		}
	}
	random_device seed;
	double A,B,C,x,y,z;
	A=seed(),B=seed(),C=seed();
	//for(auto &a : ions) {x=a.qx,y=a.qy,z=a.qz; ROTATION(a.qx,a.qy,a.qz,A,B,C,x,y,z);}
}

//------------------------------------------------------------------------
void
Variables::set_initial_velocity(Physical *pp) {
	random_device seed;
	default_random_engine engine(seed());
	double msi=28.085/6.02e23/1000.0;
	normal_distribution<> dist(0.0, sqrt(kb*T/msi));
	mt19937 mt(seed());
	for(auto &a : ions) {
		a.px=dist(engine)*1e-5;
		a.py=dist(engine)*1e-5;
		a.pz=dist(engine)*1e-5;
	}
}
