#pragma once
#include "variables.hpp"
#include "observer.hpp"
#include "PhysicalProp.hpp"
#include "tersoff.hpp"
#include "sw.hpp"
#include "flags.hpp"
#define N 50000000
#define M 1000
#define MM 100
//------------------------------------------------------------------------

class MD {
	private:


    public:
        double kappaIonDipole;
		Variables *vars;
		Observer *obs;
		Physical *pp;
		Tersoff *ters;
		SW *sw;
		FLAG *flags;

	//	vectors for pairlist
		std::vector<int> gas_in;	/*	gas list around ion1	*/
		std::vector<int> gas_out;	/*	gas list far from ion1	*/
		std::vector<Pair> pairs;	/*	pair list	*/
		std::vector<Pair> pairs_tersoff;	/*	pair list	*/
		std::vector<Pair> pairs_gi;	/*	gas-ion interaction pair list	*/
		std::vector<Pair> pairs_gg;	/*	gas-gas interaction pair list	*/
        std::vector<Pair> pairs_vi;	/*	vapor-ion interaction pair list	*/
		double margin_length;	

	//	velocity verlet
		void run_diff(char** argv);
		void run_recomb_pre(void);
		void run_recomb_repre(void);
		void verlet(void);
        void v_verlet(void);
        void v_verlet_eq(void);
		int verlet_recomb(MD *md2);	
		int verlet_recomb_after(MD *md2);	
		void update_position(void);		
        void update_position_v(void);		
        void update_position_constrained(void);		
		void velocity_calculation(void);
        void velocity_calculation_v(void);
        void velocity_calculation_veq(void);
        void update_velocity_constrained(void);


	//	pair list
		void make_pair(void);
		void make_pair_gasgas(void);
		void make_pair_gasion(void);
		void check_pairlist(void);
        void makeDiatomicProp_in(int i);
        void makeDiatomicProp_out(int i);

	//	initialization
		void initialization_gas(void);
        void initialization_vapor(void);
		void reinitialization(void);
		int initial_recomb(MD *md2, MD *Md1, MD *Md2, Physical *pp12);

	//	interaction
		void compute_inter(void);
		void compute_intra(void);
		int compute_domdom(MD *md2);
		int compute_domdom_after(MD *md2);
		void bond_harmonic(void);
		void angle_harmonic(void);
		void dihedral_fourier(void);
		void compute_lj_coul(void);
		void compute_lj(void);
		void compute_ion_dipole(void);
        void compute_ion_dipole2(void);
        void compute_ion_dipole3(void);
        void compute_ion_dipole4(void);
		void compute_coul(void);
		void compute_domdom_lj_coul(MD *md2);
		void compute_domdom_lj(MD *md2);
		void compute_gasgas_lj(void);
		void compute_cg(void);
        void computeN2(void);

	//	periodic
		void periodic(void);	/*	periodic condition for gas_in	*/
		void boundary_scaling_gas_move(int i);
		void boundary_scaling_ion_move(int i);
		int number;	/*	number of used distribution in fvlist	*/
		int number_N2;
		std::vector<double> vflux;	/*	seed of vf(v)/c	*/
		std::vector<N2_relative> relative_N2s;	/*	seed of vf(v)/c	*/
		int* dis;	/*	number of gas molecuole v(dv*i)<v<v(dv*(i+1)) */
		int loop, loop_update;	/*	current fixing time(loop) and update fixing time(loop) of out_gas for multi-timestep	*/
		void making_trans_dist(void);	/*	making distribution of vf(v)/c	*/
		void making_intra_dist(void);
		double pre_ion[3];


	//	analysis (calculating position and velocity of center of mass)
		void analysis_gas(void);	/*	calculation of center of ion1 and ion2, also collision judgement of collision and not collision	*/
		void analysis_ion(void);	/*	calculation of center of ion1 and ion2, also collision judgement of collision and not collision	*/
		double ion_r[3];
		double ion_v[3];	/*	center of ion1 and ion2	*/
		double ion_f[3];
		double gas_r[3];
		double gas_v[3];	/*	center of ion1 and ion2	*/
		double gyration;

	//	gas velocity distribution
		void velocity_distribution(void);
		void update_distribution(void);
		double disv[MM];
		double disvx_pos[MM];
		double disvy_pos[MM];
		double disvz_pos[MM];
		double disvx_neg[MM];
		double disvy_neg[MM];
		double disvz_neg[MM];
		const int UPDATE_DISTRIBUTION = 1000000000;
		const int VELOCITY_DISTRIBUTION = 100000;

	//	rigid
		int verlet_rigid(MD *md2);	
		void rigid_to_flex(MD *md2, MD *Md1, MD *Md2);
		void flex_to_rigid(MD *md2, MD *Md1, MD *Md2);
		int compute_domdom_rigid(MD *md2);
		double rigid_dt;
	
	//	coarse grain
		void verlet_cg(void);
		void velocity_calculation_cg(void);
		void update_position_cg(void);
		int verlet_recomb_cg(MD *md2);

	//	export
		void export_dump_recomb(MD *md2);
		void export_dump(void);
		void export_dump_CG(void);
		void export_dump_recomb_CG(MD *md2);
		void export_dump_rigid(MD *md2);

		/*other*/
		void output(void);
		void output_gas(void);
		void output_temp(double gastemp, double iontemp);
		void output_initial(void);
		void display(int output_ONOFF);	
		char filepath[100];
		char filepath_gyration[100];
		void fix_cell_center(void);
		void gyration_initial(void);
		void gyration_out(MD *md2);
		string gyration_path;
		double crsq;

		void velocity_scaling(void);	
		void nosehoover_ion(void);
		void nosehoover_zeta(void);
		void nosehoover_gas(void);
		void nosehoover_zeta_gas(void);
		double zeta;
		void setNVE(void);
		void setNVErecomb(void);
		void setNVTion(double temp);

		double del2,CD2,rmin2;
		double v0, th0;


		MD(char* condfile, char* file);
		~MD(void);
		void run(char** argv);
		int yesno;	/*	flag for collision or not collision	*/
};

//------------------------------------------------------------------------
