//------------------------------------------------------------------------
#include "CalculationCond.hpp"
#include "PhysicalProp.hpp"
#include "md.hpp"
#include "run.hpp"
//------------------------------------------------------------------------



int main(int argc,char *argv[]) {
/////////////////////////////////////////////////////////////////////
/*	
	Recombination simulation
	-Diffusion coefficient and collision rate coefficient simulation 
	require 2 and 7 input parameters, respectively.
	- 1:molecular infomation file name 
	- 2:calculation condition file name 
	- 3:2nd molecular infomation file name 
	- 4:2nd calculation condition file name 
	- 5:coefficient of delta 
	- 6:collision distance 
	- 7:calculation number
*/
/////////////////////////////////////////////////////////////////////
	if(argc==8){
		calculation_number=stoi(argv[7]);
//	output file initialization
		char output_file[100];
		sprintf(output_file, "delta%s_distance%s_z1_%s.dat", argv[5], argv[6], argv[7]); 
		FILE*ff=fopen(output_file, "w");
		fclose(ff);

//	construct MD class
		MD *Md1 = new MD(argv[2], argv[1]);	
		MD *Md2 = new MD(argv[4], argv[3]);
		Md1->making_trans_dist();
		Md2->making_trans_dist();
		Md1->run_recomb_pre();
		Md2->run_recomb_pre();

		FILE*f=fopen(Md1->pp->dump_path, "w");
		fclose(f);

//	set physical relative properties (reduced mass, diffusivity etc...)
//	Md classes are constructed for only pre-calculation (setting of each
//	simulation.
		Physical *pp12=new Physical();
		pp12->set_recomb(Md1->pp, Md2->pp,double(atof(argv[5])),double(atof(argv[6])));

//	recombinatin calculation loop
//	md classes are constructed fro main recombinaiton simulation.
//	after finishing a calculation, the classes are deleted and new classes  
//	are re-constructed.
		int i=0;	//	N collision
		int total=0;	//	N calculation
		while(i<100){
			ff=fopen(output_file, "a");	
			Md1->run_recomb_repre();	//	re-preparation of the simulation
			Md2->run_recomb_repre();	//	re-preparation of the simulation
			MD *md1 = new MD(argv[2], argv[1]);
			MD *md2 = new MD(argv[4], argv[3]);
			int ii=-1;
			total+=md1->initial_recomb(md2,Md1,Md2,pp12);
			ii=run_recomb(md1, md2, Md1, Md2);
			if(ii==1) {
				i++;
				run_after(md1,md2,argv[7],i);
				md1->gyration_out(md2);
			}

			fprintf(ff, "%d\t%d\t%f\t%d\t%f\t%f\n", i, total, md1->rmin2, ii, pp12->v0, pp12->th0);
			fclose(ff);
			delete md1;
			delete md2;
		}
	}



/////////////////////////////////////////////////////////////////////
/*	
	Diffusion coefficient estimation	
	-Diffusion coefficient and collision rate coefficient simulation 
	require 2 and 7 input parameters, respectively.
	- 1:molecular infomation file name 
	- 2:calculation condition file name 
	- 3:calculation number
*/
/////////////////////////////////////////////////////////////////////

	if(argc==4){
		calculation_number=stoi(argv[3]);
		MD *md = new MD(argv[2], argv[1]);	
		md->making_trans_dist();
		md->run_diff(argv);
	}


/////////////////////////////////////////////////////////////////////
/*	
	Exception error
*/
/////////////////////////////////////////////////////////////////////
	if(argc!=8&&argc!=4) cout<<"Error:Number of input parameters. -> Diffusion coefficient and collision rate coefficient simulation require 2 and 7 input parameters, respectively \n1:molecular infomation file name \n2:calculation condition file name \n3:2nd molecular infomation file name \n4:2nd calculation condition file name \n5:coefficient of delta \n6:collision distance \n7:calculation number"<<endl;
	return 0;
}

