#pragma once
#include "CalculationCond.hpp"
#include "flags.hpp"

class Physical{
public:
	const double dN2=3.1e-10;	/*	diameter	*/
	const double MN2=28.0, MHe=4.0026, MAr=39.948;
	const double mN2=MN2/Nw/1000.0;
	const double myuHe=1.47e-5, myuN2=1.7e-5, myuAr=2.28e-5;
	const double alphaHe=0.208, alphaN2=1.7, alphaAr=1.664;
	double Mion, Mgas, mgas,m , m_gas, m_N2, z, myu, D, D0, ram, delta, alphagas,v0,th0,CD,D0_He,D0_Ar,D0_N2,zzee;
	int num_gas;
	double rm, b, bm, integral, guzai, ita, itafm, itac, cgas;
	double d;	/* diameter of ions */
	double friction;
	double c;
	char *gyration_path;
	char *dump_path;
	char *RDF_path;
	double Tnh_gas, Tnh_ion;
	int ft;

	long int step_relax;
	long int step_repre;

	/*******************FUNCTION************************/
	double DELTA(double ramion, double ap);
	double RAM(double D);
	double DIFFUSION(double ap);
	double OME(double mred, double D);
	void PhysicalProp_set(char* condfile, char* file, FLAG *flags);
	void set_recomb(Physical *pp1, Physical *pp2, double Cdelta, double CD_in);
	double thc(double v);
	double CG_r0;

private:
	void read(char* infile);
	void set_condition(char* condfile, FLAG *flags);
	std::vector<double> atype;
};

