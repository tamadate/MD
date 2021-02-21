//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*	
	make pair list
	- always gas-ion pair list was calculated
	- if gas-gas interaction is ON, make pair_gasgas
	- get max velocity of gas molecule vmax
	- set update loop = margine_length/(vmax*dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::make_pair(void){
	make_pair_gasion();
	if(flags->force_gasgas_lj==1) make_pair_gasgas();
//	set number of steps to update pair list
	Gas *gases = vars->gases.data();
	double vmax2 = 0.0;
	int gos=gas_out.size();
	for (auto &a : vars->gases) {
		double px=a.px;
		double py=a.py;
		double pz=a.pz;
		double v2 = px*px + py*py + pz*pz;
		if (vmax2 < v2) vmax2 = v2;
	}
	double vion2=ion_v[0]*ion_v[0]+ion_v[1]*ion_v[1]+ion_v[2]*ion_v[2];
	if (vmax2<vion2) vmax2=vion2;
	double vmax = sqrt(vmax2)*2.0;
	loop_update=margin_length/(vmax*dt);
	loop=0;
}


/////////////////////////////////////////////////////////////////////
/*	
	make gas-ion pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::make_pair_gasion(void){
// clear gas_in, gas_out and pair list of gas-ion
	gas_in.clear();
	gas_out.clear();
	pairs_gi.clear();

	Gas *gases = vars->gases.data();
	Ion *ions = vars->ions.data();
	int is=vars->ions.size();
	int gs=vars->gases.size();
	if(flags->cg==0){
		for (int i=0;i<gs;i++){
			int flag_in=0;
			for (int j=0;j<is;j++){
				double dx = gases[i].qx - ions[j].qx;
				double dy = gases[i].qy - ions[j].qy;
				double dz = gases[i].qz - ions[j].qz;
				double r2 = (dx * dx + dy * dy + dz * dz);
				if (r2 < ML2){
					flag_in++;
					Pair p;
					p.i=i;
					p.j=j;
					pairs_gi.push_back(p);
				}
			}
	//if inter-gas interaction flag is ON, stand flag_in ON
	//if flag_in ON, molecule push back to gas_in
	//if flag_in OFF, molecule push back to gas_out
			if(flags->force_gasgas_lj==1) flag_in=1;	
			if(flag_in>0) gas_in.push_back(i);
			if(flag_in==0) gas_out.push_back(i);
		}
	}
	if(flags->cg==1){
		for (int i=0;i<gs;i++){
			int flag_in=0;
			double dx = gases[i].qx - ion_r[0];
			double dy = gases[i].qy - ion_r[1];
			double dz = gases[i].qz - ion_r[2];
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < 1600)	gas_in.push_back(i);
			else	gas_out.push_back(i);
		}
	}
}


/////////////////////////////////////////////////////////////////////
/*	
	make gas-gas pair list
*/
/////////////////////////////////////////////////////////////////////
void
MD::make_pair_gasgas(void){
	pairs_gg.clear();
	Gas *gases = vars->gases.data();
	int gs=vars->gases.size();
	for (int i=0; i<gs-1; i++){
		for (int j=i+1; j<gs; j++){
			double dx = gases[i].qx - gases[j].qx;
			double dy = gases[i].qy - gases[j].qy;
			double dz = gases[i].qz - gases[j].qz;
			adjust_periodic(dx, dy, dz);
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				pairs_gg.push_back(p);
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////
/*	
	check necessity of the pair list updating
*/
/////////////////////////////////////////////////////////////////////
void
MD::check_pairlist(void){
	loop++;
	if(loop>loop_update){
		periodic_out(0,gas_out);
		int gos=gas_out.size();
		Gas *gases = vars->gases.data();
		for (int k = 0; k < gos; k++) {
			double dqx=0,dqy=0,dqz=0;
			for(int i = gas_out[k]; i < gas_out[k]+pp->num_gas; i++){
				dqx+=gases[i].px*gases[i].mass;
				dqy+=gases[i].py*gases[i].mass;
				dqz+=gases[i].pz*gases[i].mass;
			}
			dqx=dqx*dt*loop/pp->Mgas;
			dqy=dqy*dt*loop/pp->Mgas;
			dqz=dqz*dt*loop/pp->Mgas;
			for(int i = gas_out[k]; i < gas_out[k]+pp->num_gas; i++){
				gases[i].qx += dqx;
				gases[i].qy += dqy;
		  		gases[i].qz += dqz;
			}
		}
		periodic_out(1,gas_out);
		make_pair();
		for(int i=0;i<3;i++) pre_ion[i]=ion_r[i];
	}
	if(flags->force_sw==1) sw->check_pairlist(vars);
	if(flags->force_ters==1) ters->check_pairlist(vars);
}

