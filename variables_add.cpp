//------------------------------------------------------------------------
#include "variables.hpp"
//------------------------------------------------------------------------
void
Variables::add_atype(int type, double mass, double coeff1, double coeff2) {
 	Atom_type at;
	at.type = type;
	at.mass = mass; 
	at.coeff1 = coeff1;
	at.coeff2 = coeff2;
	atypes.push_back(at);
}
void
Variables::add_btype(int type, double coeff1, double coeff2) {
 	Bond_type bt;
	bt.type = type;
	bt.coeff1 = coeff1;
	bt.coeff2 = coeff2;
	btypes.push_back(bt);
}
void
Variables::add_ctype(int type, double coeff1, double coeff2) {
 	Angle_type ct;
	ct.type = type;
	ct.coeff1 = coeff1;
	ct.coeff2 = coeff2;
	ctypes.push_back(ct);
}
void
Variables::add_dtype(int type, int multi, double coeff1, double coeff2, double coeff3, double coeff4, double coeff5, double coeff6, double coeff7, double coeff8, double coeff9, double coeff10, double coeff11, double coeff12, double coeff13, double coeff14, double coeff15) {
 	Dihedral_type dit;
	dit.type = type;
	dit.coeff1 = coeff1;
	dit.coeff2 = coeff2;
	dit.coeff3 = coeff3;
	dit.coeff4 = coeff4;
	dit.coeff5 = coeff5;
    dit.coeff6 = coeff6;
    dit.coeff7 = coeff7;
    dit.coeff8 = coeff8;
    dit.coeff9 = coeff9;
    dit.coeff10 = coeff10;
    dit.coeff11 = coeff11;
    dit.coeff12 = coeff12;
    dit.coeff13 = coeff13;
    dit.coeff14 = coeff14;
    dit.coeff15 = coeff15;
	dtypes.push_back(dit);
}
void
Variables::add_ions(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass) {
 	Ion a;
	a.id = id;
	a.type = type; 
	a.qx = x;
	a.qy = y;
	a.qz = z;
	a.px = vx;
	a.py = vy;
	a.pz = vz;
	a.mass = mass;
	a.fx = fx;
	a.fy = fy;
	a.fz = fz;
	a.charge = charge;
	a.ix=0;
	a.iy=0;
	a.iz=0;
	ions.push_back(a);
}
void
Variables::add_gases(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass) {
 	Gas a;
	a.id = id;
	a.type = type; 
	a.qx = x;
	a.qy = y;
	a.qz = z;
	a.px = vx;
	a.py = vy;
	a.pz = vz;
	a.mass = mass;
	a.fx = fx;
	a.fy = fy;
	a.fz = fz;
	a.charge = charge;
	a.ix=0;
	a.iy=0;
	a.iz=0;
	gases.push_back(a);
}
void
Variables::add_vapors(int id, int type, double x, double y, double z, double vx, double vy, double vz, double fx, double fy, double fz, double charge, double mass) {
 	Atom a;
	a.id = id;
	a.type = type; 
	a.qx = x;
	a.qy = y;
	a.qz = z;
	a.px = vx;
	a.py = vy;
	a.pz = vz;
	a.mass = mass;
	a.fx = fx;
	a.fy = fy;
	a.fz = fz;
	a.charge = charge;
	a.ix=0;
	a.iy=0;
	a.iz=0;
	vapors.push_back(a);
}
void
Variables::add_bonds(int atom1, int atom2, int type) {
 	Bond b;
	b.atom1 = atom1;
	b.atom2 = atom2; 
	b.type = type;
	bonds.push_back(b);
}
void
Variables::add_angles(int atom1, int atom2, int atom3, int type) {
 	Angle c;
	c.atom1 = atom1;
	c.atom2 = atom2; 
	c.atom3 = atom3;
	c.type = type;
	angles.push_back(c);
}
void
Variables::add_dihedrals(int atom1, int atom2, int atom3, int atom4, int type) {
 	Dihedral d;
	d.atom1 = atom1;
	d.atom2 = atom2;
	d.atom3 = atom3;
	d.atom4 = atom4;
	d.type = type;
	dihedrals.push_back(d);
}
