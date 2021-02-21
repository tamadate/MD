//------------------------------------------------------------------------
#include "PhysicalProp.hpp"
//------------------------------------------------------------------------


#define SIZE 1000000

int gastype=2;	/*1:He, 2:Ar, 3:N2*/
long int Noftimestep=10000000000;
double p=pow(10,(2.0+3.0));	/*Pa*/
double T = 300;     /* temp. */
int calculation_number=7;
int Nof_around_gas=1000;
double d_size=pow(Nof_around_gas*kb*T/p,1/3.0)*1e10;
double V=d_size*d_size*d_size;
double dt = 0.5;	/*	fs	*/
double CUTOFF = 10.0;	/*	A	*/
double MARGIN = 15.0;	/*	A	*/
double ML2 = (CUTOFF+MARGIN)*(CUTOFF+MARGIN);
double CL2 = (CUTOFF*CUTOFF);
double MARGIN_PEG = 3.0;	/*	A	*/
double ML2_PEG = (CUTOFF+MARGIN)*(CUTOFF+MARGIN);

void
Physical::PhysicalProp_set(char* condfile, char* file, FLAG *flags){	
	set_condition(condfile, flags);
	read(file);
	if(gastype==1) {
		Mgas=MHe;
		D0=D0_He;
		myu=myuHe;
		alphagas=alphaHe;
		num_gas=1;
	}
	if(gastype==2) {
		Mgas=MN2;
		D0=D0_N2;
		myu=myuN2;
		alphagas=alphaN2/2.0;
		num_gas=2;
	}
	if(gastype==3) {
		Mgas=MN2;
		D0=D0_N2;
		myu=myuN2;
		alphagas=alphaN2;
		num_gas=1;
	}
	if(gastype==4) {
		Mgas=MAr;
		D0=D0_Ar;
		myu=myuAr;
		alphagas=alphaAr;
		num_gas=1;
	}
	mgas=Mgas/Nw/1000.0;
	m=Mion/Nw/1000.0;
	m_gas=mgas*m/(mgas+m);
	m_N2=mN2*m/(mN2+m);
	d=pow(OME(m_N2,D0_N2)*4.0/1.36/M_PI,0.5)-dN2;
	D=D0*1e5/p;//DIFFUSION((pow(OME(m_gas,D0)*4.0/1.36/M_PI,0.5))/2.0);
	ram=RAM(D);
	friction=kb*T/D;
	delta = DELTA(ram,d/2.0);
	c=sqrt(8*kb*T/M_PI/m);
	cgas=sqrt(8*kb*T/M_PI/mgas);

	printf("Delta\t\t\t%f ang.\nLambda\t\t\t%f ang.\n",delta*1e10,ram*1e10);
	printf("Ion mass\t\t%f g/mol\nIon charges\t\t%f\nIon diameter\t\t%f ang.\n", Mion,z,d*1e10);
	cout<<"**************************************************"<<endl;
	cout<<"**************************************************"<<endl;
	cout<<"**************************************************"<<endl;
}


void
Physical::read(char* infile){
	ifstream stream(infile);
	string str;
	int iflag=0;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="atom type name mass coeff1 coeff2") {iflag=1; continue;}
		if (str=="atoms") {iflag=2; continue;}
		if (str=="bonds") {break;}
	    string tmp;
    	istringstream stream(str);		
		if (iflag==1) {
			int loop=0;
			double mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==2) mass=stod(tmp);
				loop++;
			}
			atype.push_back(mass);
		}
		if (iflag==2) {
			int loop=0;
			int type, id;
			double charge, x, y, mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==1) type=stoi(tmp);
				if (loop==2) charge=stod(tmp);
				loop++;
			}
			mass=atype[type-1];
			Mion+=mass;
			z+=charge;
		}
	}	
}


void
Physical::set_condition(char* condfile, FLAG *flags){
	ifstream stream(condfile);
	string str;
	int iflag=0;
	cout<<"**************************************************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"************  Calculation conditions  ************"<<endl;
	cout<<"************--------------------------************"<<endl;
	cout<<"**************************************************"<<endl;
	while(getline(stream,str)) {
		if(str.length()==0) continue;
		if (str=="Gas type") {iflag=1; continue;}
		if (str=="Calculation number") {iflag=2; continue;}
		if (str=="Temperature"||str=="T") {iflag=3; continue;}
		if (str=="Pressure"||str=="P") {iflag=4; continue;}
		if (str=="Time step"||str=="dt") {iflag=5; continue;}
		if (str=="Total number of steps") {iflag=6; continue;}
		if (str=="Number of steps for relax") {iflag=7; continue;}
		if (str=="Number of around gas") {iflag=8; continue;}
		if (str=="Cut off length") {iflag=9; continue;}
		if (str=="Margin size") {iflag=10; continue;}	
		if (str=="Diffusion coefficient in He") {iflag=11; continue;}	
		if (str=="Diffusion coefficient in N2") {iflag=12; continue;}	
		if (str=="Number of steps for repre") {iflag=13; continue;}	
		if (str=="Collision distance") {iflag=14; continue;}	
		if (str=="Diffusion coefficient in Ar") {iflag=15; continue;}
		if (str=="Gyration path") {iflag=16; continue;}	
		if (str=="RDF path") {iflag=17; continue;}	
		if (str=="gas gas interaction") {iflag=18; continue;}	
		if (str=="gas ion interaction") {iflag=19; continue;}	
		if (str=="ion ion interaction") {iflag=20; continue;}	
		if (str=="Nose-Hoover ion") {iflag=21; continue;}	
		if (str=="Nose-Hoover gas") {iflag=22; continue;}	
		if (str=="Output") {iflag=23; continue;}	
		if (str=="Coarse grain") {iflag=24; continue;}	
		if (iflag==1) {
			gastype=stoi(str);
			if(stoi(str)==1) cout<<"Gastype\tHe"<<endl;
			if(stoi(str)==2) cout<<"Gastype\tN2(di)"<<endl;
			if(stoi(str)==3) cout<<"Gastype\tN2(mono)"<<endl;
			if(stoi(str)==4) cout<<"Gastype\tAr"<<endl;
		}
		if (iflag==2) {
			calculation_number=stoi(str);
		}
		if (iflag==3) {
			T=stod(str);
			cout<<"Temperature\t\t"<<T<<" K"<<endl;
		}
		if (iflag==4) {
			p=stod(str);
			cout<<"Pressure\t\t"<<p<<" Pa"<<endl;
		}
		if (iflag==5) {
			dt=stod(str);
			cout<<"Time step\t\t"<<dt<<" fs"<<endl;
		}
		if (iflag==6) {
			double inter=stod(str);
			Noftimestep=inter;
			cout<<"Total steps\t\t"<<float(Noftimestep)<<endl;
		}
		if (iflag==7) {
			double inter=stod(str);
			step_relax=inter;
			cout<<"Relax steps\t\t"<<float(step_relax)<<endl;
		}
		if (iflag==8) {
			Nof_around_gas=stoi(str);
			cout<<"Number of gases\t\t"<<Nof_around_gas<<endl;
		}
		if (iflag==9) {
			CUTOFF=stod(str);
			cout<<"Cutoff\t\t\t"<<CUTOFF<<" ang."<<endl;
		}
		if (iflag==10) {
			MARGIN=stod(str);
			cout<<"Margin size\t\t"<<MARGIN<<" ang."<<endl;
		}
		if (iflag==11) {
			D0_He=stod(str);
			cout<<"Diffusion coeff. in He\t"<<D0_He<<" cm^2 s^-1."<<endl;
		}
		if (iflag==12) {
			D0_N2=stod(str);
			cout<<"Diffusion coeff. in N2\t"<<D0_N2<<" cm^2 s^-1."<<endl;
		}
		if (iflag==13) {
			step_repre=stoi(str);
			cout<<"Re-relax steps\t\t"<<float(step_repre)<<endl;
		}
		if (iflag==14) {
			/*CD=stod(str);*/
		}
		if (iflag==15) {
			D0_Ar=stod(str);
			cout<<"Diffusion coeff. in Ar\t"<<D0_Ar<<" cm^2 s^-1."<<endl;
		}
		if (iflag==16) {
			ostringstream ss;
			ss<<str<<"_"<<calculation_number<<".dat";
			string tmp=ss.str();
			gyration_path=new char[tmp.length()+10];
			strcpy(gyration_path,tmp.c_str());
			flags->gyration=1;
			FILE*f=fopen(gyration_path, "w");
			fclose(f);
			cout<<"Gyration --> ON -->\t"<<tmp<<endl;
		}
		if (iflag==17) {
			RDF_path=new char[str.length()+1];
			strcpy(RDF_path,str.c_str());
			flags->RDF=1;
			cout<<"RDF\tON\tPath="<<str<<endl;
		}
		if (iflag==18) {
			if (str=="LJ") {
				flags->force_gasgas_lj=1;
				cout<<"gas-gas interaction -->\tON"<<endl;
			}
			else if (str=="LJsemiNVT") {
				flags->force_gasgas_lj=1;
				flags->semi_NVT_gasgas=1;
				cout<<"gas-gas interaction -->\tON"<<endl;
				cout<<"Boundary rescaling -->\tON"<<endl;
			}
			else if (str=="OFF") {
				flags->force_gasgas_lj=0;
				flags->semi_NVT_gasgas=0;
				cout<<"gas-gas interaction -->\tOFF"<<endl;
				cout<<"Boundary rescaling -->\tOFF"<<endl;
			}
			else printf("**************Uknown gas gas parameter was found**************\n");
		}
		if (iflag==19) {
			if (str=="LJ") flags->force_lj=1;
			else if (str=="ion dipole") flags->force_ion_dipole=1;
			else if (str=="OFF") {
				flags->force_ion_dipole=0;
				flags->force_lj=0;
			}
			else printf("**************Uknown gas ion parameter was found**************\n");
		}
		if (iflag==20) {
			if (str=="AMBER") flags->intra_AMBER=1;
			else if (str=="coul") flags->force_coul=1;
			else if (str=="Stilinger-Weber") flags->force_sw=1;
			else if (str=="Tersoff") flags->force_ters=1;
			else printf("**************Uknown ion ion parameter was found**************\n");
		}
		if (iflag==21) {
			if (str=="OFF") {flags->nose_hoover_ion=0; cout<<"Nose-Hoover for ion --> OFF"<<endl;}
			else {
				flags->nose_hoover_ion=1;
				Tnh_ion=stod(str);
				cout<<"Nose-Hoover for ion --> ON --> "<<Tnh_ion<<" K"<<endl;
			}
		}
		if (iflag==22) {
			if (str=="OFF") {flags->nose_hoover_gas=0; cout<<"Nose-Hoover for gas --> OFF"<<endl;}
			else {
				flags->nose_hoover_gas=1;
				Tnh_gas=stod(str);
				cout<<"Nose-Hoover for gas --> ON --> "<<Tnh_gas<<" K"<<endl;
			}
		}
		if (iflag==23) {
			ostringstream ss;
			ss<<str<<"_"<<calculation_number<<".dump";
			string tmp=ss.str();
			dump_path=new char[tmp.length()+1];
			strcpy(dump_path,tmp.c_str());
			cout<<"Dump file -->\t\t"<<tmp<<endl;
		}
		if (iflag==24) {
			istringstream stream(str);	
			string tmp;
			int loop=0;
			double mass;
			while(getline(stream,tmp,'\t')) {
				if (loop==0) flags->cg=stoi(tmp);
				if (loop==1) CG_r0=stod(tmp);
				loop++;
			}
		}
	}
	d_size=pow(Nof_around_gas*kb*T/p,1/3.0)*1e10;//pow(28.0855*8/6.02e23/(2.218e-24),1/3.0)*5;
	V=d_size*d_size*d_size;
	ML2 = (CUTOFF+MARGIN)*(CUTOFF+MARGIN);
	cout<<"Domain size\t\t"<<d_size<<" ang."<<endl;
}


void
Physical::set_recomb(Physical *pp1, Physical *pp2, double Cdelta, double CD_in){	
	Mgas=pp1->Mgas;
	myu=pp1->myu;
	alphagas=pp1->alphagas;
	num_gas=pp1->num_gas;
	mgas=pp1->mgas;
	m=1/(1/pp1->m+1/pp2->m);
	Mion=1/(1/pp1->Mion+1/pp2->Mion);
	m_gas=mgas*m/(mgas+m);
	m_N2=mN2*m/(mN2+m);
	d=(pp1->d+pp2->d)*0.5;
	D=pp1->D+pp2->D;
	ram=RAM(D);
	friction=kb*T/D;
	delta = DELTA(ram,d/2.0)*Cdelta;
	c=sqrt(8*kb*T/M_PI/m);
	cgas=sqrt(8*kb*T/M_PI/mgas);
	CD=CD_in;
	pp1->delta=pp2->delta=delta;
	pp1->CD=pp2->CD=CD;
	zzee=pp1->z*pp2->z*e*e;

	cout << "\n\n===CALCULATION CONDITION====" << endl;
	printf("TEMP=%1.2fK\nPRESS=%1.2fkPa\nDELTA=%1.2fA\nRAMDA=%1.2fA\n", T,p/1000.0,delta*1e10,ram*1e10);
	cout << "============GAS=============" << endl;
	printf("MASS=%1.2fg/mol\nVISCOSITY=%1.2ePas\n", Mgas,myu);
	cout << "============ION=============" << endl;
	printf("REDUCED MASS=%1.2fg/mol\nCHARGE=%1.2f\tvs.\t%1.2f\nDIAMETER=%1.2fA\nDIFF COEFF=%1.2fcm2/s\n", Mion,pp1->z,pp2->z,d*1e10,D*1e4);
	cout << "============================" << endl;
}


//------------------------------------------------------------------------
double 
Physical::OME(double mred, double D){
	return 3.0*kb*T/16.0*pow(2*M_PI*kb*T/mred,0.5)/1.013e5/D;
}
//------------------------------------------------------------------------
double 
Physical::RAM(double D){
	return 1.329*D*sqrt(m_gas/kb/T);
}
//------------------------------------------------------------------------
double 
Physical::DIFFUSION(double ap){
	double ramda=myu*kb*T/0.491/pow(8*kb*T/M_PI/mgas,0.5)/p/mgas;
	double Cc=1+ramda/ap*(1.257+0.4*exp(-1.1*ap/ramda));
	return kb*T/6/M_PI/myu/ap*Cc;
}
//------------------------------------------------------------------------
double 
Physical::DELTA(double ramion, double ap){
	return (pow(ap, 3.0) / pow(ramion, 2.0)) * (1 / 5.0 * pow(1 + ramion / ap, 5.0) - 1 / 3.0 * (1 + ramion * ramion / ap / ap) * pow(1 + ramion / ap, 3.0) + 2 / 15.0 * pow(1 + ramion * ramion / ap / ap, 5 / 2.0));
}
//------------------------------------------------------------------------
double 
Physical::thc(double v){
	double in_root=1-zzee*0.5/M_PI/m/v/v/e0*(1e10/CD-1/delta);
	double bc=CD*sqrt(in_root)*1e-10;
	return asin(bc/delta)/M_PI*180;
}
//------------------------------------------------------------------------


