#pragma once
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "CalculationCond.hpp"
#include "PhysicalProp.hpp"

using namespace std;
//------------------------------------------------------------------------
struct Atom_type {
 	int type;
	double mass;
	double coeff1, coeff2;
};
//------------------------------------------------------------------------
struct Bond_type {
 	int type;
	double coeff1, coeff2;
};
//------------------------------------------------------------------------
struct Angle_type {
 	int type;
	double coeff1, coeff2;
};
//------------------------------------------------------------------------
struct Dihedral_type {
 	int type, multi;
	double coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8, coeff9, coeff10,  coeff11, coeff12, coeff13, coeff14, coeff15;
};
//------------------------------------------------------------------------
struct Ion {
 	double qx, qy, qz;
	double px, py, pz;
	double charge;
	int type;
	int id;
	double fx, fy, fz;
	double mass;
	int ix, iy, iz;
};
//------------------------------------------------------------------------
struct Gas {
 	double qx, qy, qz;
	double px, py, pz;
	double charge;
	int type;
	int id;
	double fx, fy, fz;
	double mass;
	int ix, iy, iz;
};
//------------------------------------------------------------------------
struct Atom {
 	double qx, qy, qz;
	double px, py, pz;
	double charge;
	int type;
	int id;
	double fx, fy, fz;
	double mass;
	int ix, iy, iz;
};
//------------------------------------------------------------------------
struct Bond {
	int atom1, atom2, type;
};
//------------------------------------------------------------------------
struct Angle {
	int atom1, atom2, atom3, type;
};
//------------------------------------------------------------------------
struct Dihedral {
	int atom1, atom2, atom3, atom4, type;
};
//------------------------------------------------------------------------
struct Pair{
	int i,j;
};
struct Pair_many{
	int i;
	std::vector<int> j;	/*	pair list	*/
};

//------------------------------------------------------------------------
struct N2_relative{
	double x,y,z,vx1,vy1,vz1,vx2,vy2,vz2;
};
//------------------------------------------------------------------------
class Variables {
public:

	Variables(void) {time = 0.0;}

	/*variables*/
	std::vector<Ion> ions;
	std::vector<Gas> gases;
    std::vector<Atom> vapors;
	double time;
	double zeta_ion;
	double zeta_gas;

	/*vectors for potential calculation*/
    std::vector<Pair> ion_pairs;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;
	std::vector<Atom_type> atypes;
	std::vector<Bond_type> btypes;
	std::vector<Angle_type> ctypes;
	std::vector<Dihedral_type> dtypes;
	std::vector<int> molecules;
	std::vector<vector<vector<double>>> pair_coeff;

	/*add to vectors*/
	void add_ions(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass);
	void add_gases(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass);
    void add_vapors(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass);
	void add_bonds(int atom1, int atom2, int type);
	void add_angles(int atom1, int atom2, int atom3, int type);
	void add_dihedrals(int atom1, int atom2, int atom3, int atom4, int type);
	void add_atype(int type, double mass, double coeff1, double coeff2);
	void add_btype(int type, double coeff1, double coeff2);
	void add_ctype(int type, double coeff1, double coeff2);
	void add_dtype(int type, int multi, double coeff1, double coeff2, double coeff3, double coeff4, double coeff5, double coeff6, double coeff7, double coeff8, double coeff9, double coeff10, double coeff11, double coeff12, double coeff13, double coeff14, double coeff15);

	/*initialization and export to dump file*/
	void read_initial(char* infile);
	void set_initial_velocity(Physical *pp);

	void ROTATION(double &X, double &Y, double &Z, double A, double B, double C, double x, double y, double z);
private:
};
//------------------------------------------------------------------------
