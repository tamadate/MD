#include "md.hpp"
#include <random>
#include <algorithm>
#include <dirent.h>
#include <vector>
/*########################################################################################

-----Periodic conditions-----

#######################################################################################*/

/**************************ajust periodic**********************************/
void
adjust_periodic(double &dx, double &dy, double &dz) {
	const double LH = d_size * 0.5;
	if (dx < -LH)dx += d_size;
	if (dx > LH) dx -= d_size;
	if (dy < -LH)dy += d_size;
	if (dy > LH) dy -= d_size;
	if (dz < -LH)dz += d_size;
	if (dz > LH) dz -= d_size;
}
/************************periodec condition for out gas******************************/
void
MD::periodic_out(int dist, std::vector<int> &gas) {
	Gas *gases = vars->gases.data();
	int gos=gas.size();
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	for (int k=0;k<gos;k++) {
		int i=gas[k];
		flag=flagx=flagy=flagz=0;
		if (gases[i].qx < ion_r[0]-HL) gases[i].qx += d_size, flagx--, flag++;
		if (gases[i].qy < ion_r[1]-HL) gases[i].qy += d_size, flagy--, flag++;
		if (gases[i].qz < ion_r[2]-HL) gases[i].qz += d_size, flagz--, flag++;
		if (gases[i].qx > ion_r[0]+HL) gases[i].qx -= d_size, flagx++, flag++;
		if (gases[i].qy > ion_r[1]+HL) gases[i].qy -= d_size, flagy++, flag++;
		if (gases[i].qz > ion_r[2]+HL) gases[i].qz -= d_size, flagz++, flag++;
		if (i%pp->num_gas==0&&flag>0&&dist==1) {
			if(number>vflux.size()*0.8) {making_trans_dist();}
			if(number_N2>relative_N2s.size()) {making_intra_dist();}
			vx=0,vy=0,vz=0;
			for(int l = i; l < i+pp->num_gas; l++){
				vx+=gases[l].px*gases[l].mass;
				vy+=gases[l].py*gases[l].mass;
				vz+=gases[l].pz*gases[l].mass;
			}
			vx = vx/pp->m_gas;
			vy = vy/pp->m_gas;
			vz = vz/pp->m_gas;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = vflux[number]*1e-5;
			mod_factor= vMB/v;
			if(gastype==2){
				gases[i].px = vx*mod_factor+(relative_N2s[number_N2].vx1-relative_N2s[number_N2].vx2)*0.5;
				gases[i].py = vy*mod_factor+(relative_N2s[number_N2].vy1-relative_N2s[number_N2].vy2)*0.5;
				gases[i].pz = vz*mod_factor+(relative_N2s[number_N2].vz1-relative_N2s[number_N2].vz2)*0.5;
				gases[i+1].px = vx*mod_factor-(relative_N2s[number_N2].vx1-relative_N2s[number_N2].vx2)*0.5;
				gases[i+1].py = vy*mod_factor-(relative_N2s[number_N2].vy1-relative_N2s[number_N2].vy2)*0.5;
				gases[i+1].pz = vz*mod_factor-(relative_N2s[number_N2].vz1-relative_N2s[number_N2].vz2)*0.5;
				gases[i+1].qx = gases[i].qx+relative_N2s[number_N2].x;
				gases[i+1].qy = gases[i].qy+relative_N2s[number_N2].y;
				gases[i+1].qz = gases[i].qz+relative_N2s[number_N2].z;
				number++;
				number_N2++;
			}
			else{
				vx = gases[i].px;
				vy = gases[i].py;
				vz = gases[i].pz;
				v2 = vx*vx+vy*vy+vz*vz;
				v = sqrt(v2);
				vMB = vflux[number]*1e-5;
				mod_factor= vMB/v;
				gases[i].px = vx * mod_factor;
				gases[i].py = vy * mod_factor;
				gases[i].pz = vz * mod_factor;
				number++;
			}
		}
		if (i%pp->num_gas==0&&flag>0&&dist==0) {
			gases[i].px = distgas(engine) *1e-5;
			gases[i].py = distgas(engine) *1e-5;
			gases[i].pz = distgas(engine) *1e-5;
		}
	}
}
/************************periodec condition for in gas***************************/
void
MD::periodic(void) {
	Gas *gases = vars->gases.data();
	int gis=gas_in.size();
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;
	for (int k=0;k<gis;k++) {
		int i=gas_in[k];
		flag=flagx=flagy=flagz=0;
		if (gases[i].qx < ion_r[0]-HL) gases[i].qx += d_size, flagx--, flag++;
		if (gases[i].qy < ion_r[1]-HL) gases[i].qy += d_size, flagy--, flag++;
		if (gases[i].qz < ion_r[2]-HL) gases[i].qz += d_size, flagz--, flag++;
		if (gases[i].qx > ion_r[0]+HL) gases[i].qx -= d_size, flagx++, flag++;
		if (gases[i].qy > ion_r[1]+HL) gases[i].qy -= d_size, flagy++, flag++;
		if (gases[i].qz > ion_r[2]+HL) gases[i].qz -= d_size, flagz++, flag++;


	}
}

/************************periodec condition for in gas***************************/
void
MD::periodic_semi(void) {
	Gas *gases = vars->gases.data();
	int gis=gas_in.size();
	random_device seed;
	default_random_engine engine(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	double vMB, vxMB, vyMB, vzMB, vx, vy, vz, v2, v, mod_factor;

//Step 1
	for (int k=0;k<gis;k++) {
		int i=gas_in[k];
		flag=flagx=flagy=flagz=0;
		if (gases[i].qx < ion_r[0]-HL) gases[i].qx += d_size, flagx--, flag++;
		if (gases[i].qy < ion_r[1]-HL) gases[i].qy += d_size, flagy--, flag++;
		if (gases[i].qz < ion_r[2]-HL) gases[i].qz += d_size, flagz--, flag++;
		if (gases[i].qx > ion_r[0]+HL) gases[i].qx -= d_size, flagx++, flag++;
		if (gases[i].qy > ion_r[1]+HL) gases[i].qy -= d_size, flagy++, flag++;
		if (gases[i].qz > ion_r[2]+HL) gases[i].qz -= d_size, flagz++, flag++;
		if (i%pp->num_gas==0&&flag>0) {
			if(number>vflux.size()*0.8) {making_trans_dist();}
			if(number_N2>relative_N2s.size()) {making_intra_dist();}
			vx=0,vy=0,vz=0;
			for(int l = i; l < i+pp->num_gas; l++){
				vx+=gases[l].px*gases[l].mass;
				vy+=gases[l].py*gases[l].mass;
				vz+=gases[l].pz*gases[l].mass;
			}
			vx = vx/pp->m_gas;
			vy = vy/pp->m_gas;
			vz = vz/pp->m_gas;
			v2 = vx*vx+vy*vy+vz*vz;
			v = sqrt(v2);
			vMB = vflux[number]*1e-5;
			mod_factor= vMB/v;
			if(gastype==2){
				gases[i].px = vx*mod_factor+(relative_N2s[number_N2].vx1-relative_N2s[number_N2].vx2)*0.5;
				gases[i].py = vy*mod_factor+(relative_N2s[number_N2].vy1-relative_N2s[number_N2].vy2)*0.5;
				gases[i].pz = vz*mod_factor+(relative_N2s[number_N2].vz1-relative_N2s[number_N2].vz2)*0.5;
				gases[i+1].px = vx*mod_factor-(relative_N2s[number_N2].vx1-relative_N2s[number_N2].vx2)*0.5;
				gases[i+1].py = vy*mod_factor-(relative_N2s[number_N2].vy1-relative_N2s[number_N2].vy2)*0.5;
				gases[i+1].pz = vz*mod_factor-(relative_N2s[number_N2].vz1-relative_N2s[number_N2].vz2)*0.5;
				gases[i+1].qx = gases[i].qx+relative_N2s[number_N2].x;
				gases[i+1].qy = gases[i].qy+relative_N2s[number_N2].y;
				gases[i+1].qz = gases[i].qz+relative_N2s[number_N2].z;
				number++;
				number_N2++;
			}
			else{
				vx = gases[i].px;
				vy = gases[i].py;
				vz = gases[i].pz;
				v2 = vx*vx+vy*vy+vz*vz;
				v = sqrt(v2);
				vMB = vflux[number]*1e-5;
				mod_factor= vMB/v;
				gases[i].px = vx * mod_factor;
				gases[i].py = vy * mod_factor;
				gases[i].pz = vz * mod_factor;
				number++;
			}
		}
	}


	analysis_ion();

//Step 2
	for (int k=0;k<gis;k++) {
		int i=gas_in[k];
		flag=flagx=flagy=flagz=0;
		if (gases[i].qx < ion_r[0]-HL) gases[i].qx += d_size, flagx--, flag++;
		if (gases[i].qy < ion_r[1]-HL) gases[i].qy += d_size, flagy--, flag++;
		if (gases[i].qz < ion_r[2]-HL) gases[i].qz += d_size, flagz--, flag++;
		if (gases[i].qx > ion_r[0]+HL) gases[i].qx -= d_size, flagx++, flag++;
		if (gases[i].qy > ion_r[1]+HL) gases[i].qy -= d_size, flagy++, flag++;
		if (gases[i].qz > ion_r[2]+HL) gases[i].qz -= d_size, flagz++, flag++;
		if (i%pp->num_gas==0&&flag>0) {
			gases[i].px = distgas(engine) *1e-5;
			gases[i].py = distgas(engine) *1e-5;
			gases[i].pz = distgas(engine) *1e-5;
		}
	}

}

/************************making weighted MB distribution*****************************/
void
MD::making_trans_dist(void){
	number=0;
	const double vmax=5.0*pp->cgas;	/*	max velocity of gas molecule in flux velocity distribution */
	const double dv=vmax/M;	/*	resolution of flux velocity distribution */
	dis = new int[M];
	vector<double>().swap(vflux);
	for(int i=0; i<M; i++) dis[i]=0;	
	random_device seed;
	mt19937 mt(seed());
	uniform_real_distribution<double> v(0,vmax);
	int j;
	for(int i=0; i<N; i++) {
		double vnow=v(mt);
		j=vnow/dv;
		double x=(j+0.5)*dv;
		double A=pp->mgas*pp->mgas/2.0/kb/kb/T/T*x*x*x*exp(-1*pp->mgas/2.0/kb/T*x*x);
		if(dis[j]<A*30000000){
			vflux.push_back(vnow);
			dis[j]++;
		}
	}
	random_shuffle(vflux.begin(), vflux.end());
	delete [] dis;
}

void
MD::making_intra_dist(void){
	number_N2=0;
	if(gastype==2){
		vector<N2_relative>().swap(relative_N2s);
		random_device seed;
		mt19937 mt(seed());
		uniform_real_distribution<double> num(0,1000);
/*		sprintf(filepath, "N2/N2_dist_300_%d.dat",int(num(mt))); */
		sprintf(filepath, "N2_dist_300_10.dat"); 
		ifstream stream(filepath);
		string str;
		while(getline(stream,str)) {
			string tmp;
			istringstream stream(str);
			int i=0;
			N2_relative a;
			while(getline(stream,tmp,'\t')) {
				if (i==0) a.x=stod(tmp);
				if (i==1) a.y=stod(tmp);
				if (i==2) a.z=stod(tmp);
				if (i==3) a.vx1=stod(tmp);
				if (i==4) a.vy1=stod(tmp);
				if (i==5) a.vz1=stod(tmp);
				if (i==6) a.vx2=stod(tmp);
				if (i==7) a.vy2=stod(tmp);
				if (i==8) a.vz2=stod(tmp);
				i++;
			}
			double HL=d_size*0.5;
			if (a.x>HL) a.x-=d_size;
			if (a.y>HL) a.y-=d_size;
			if (a.z>HL) a.z-=d_size;
			if (a.x<-HL) a.x+=d_size;
			if (a.y<-HL) a.y+=d_size;
			if (a.z<-HL) a.z+=d_size;
			relative_N2s.push_back(a);
		}
		random_shuffle(relative_N2s.begin(), relative_N2s.end());
	}
}

