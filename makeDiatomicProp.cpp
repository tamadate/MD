//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


void 
MD::makeDiatomicProp_in(int i){
    
    if(vars->gases[i+Nof_around_gas].mass==0){
        random_device seed;
        mt19937 mt(seed());
        default_random_engine engine(seed());
        normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas*2));
        uniform_real_distribution<double> r(0,1);
        
        const double NNbond=1.098;
        const double NNbond_half=NNbond*0.5;
        //const double I1=0.5*pp->mgas*NNbond_half*NNbond_half;
        double a=r(mt)*2*M_PI;
        double b=r(mt)*2*M_PI;
        //double cosb=1-2*r(mt);
        //double b=acos(cosb);
        //double c=r(mt)*2*M_PI;
        //double wsq=-2*kb*T/I1/log(1-r(mt));
        //double w2=sqrt(wsq/(1+tan(c)*tan(c)))*1e-5;
        //double w1=tan(c)*w2;
        //cout<<sqrt(wsq)<<" "<<w1<<" "<<w2<<" "<<endl;
        
        double x1=-NNbond_half;
        double y1=0;
        double z1=0;
        double x2=NNbond_half;
        double y2=0;
        double z2=0;
        
        double vx1=0;
        double vy1=distgas(engine)*1e-5;
        double vz1=distgas(engine)*1e-5;
        double vx2=0;
        double vy2=-vy1;
        double vz2=-vz1;
        
        vars->ROTATION(x1,y1,z1,0,a,b,x1,y1,z1);
        vars->ROTATION(x2,y2,z2,0,a,b,x2,y2,z2);
        vars->ROTATION(vx1,vy1,vz1,0,a,b,vx1,vy1,vz1);
        vars->ROTATION(vx2,vy2,vz2,0,a,b,vx2,vy2,vz2);

        
        
        //  set positions
        //  index i have a center of gas molecule and i+Nof... have no information before this operation
        vars->gases[i+Nof_around_gas].qx=vars->gases[i].qx+x2;
        vars->gases[i+Nof_around_gas].qy=vars->gases[i].qy+y2;
        vars->gases[i+Nof_around_gas].qz=vars->gases[i].qz+z2;
        vars->gases[i].qx+=x1;
        vars->gases[i].qy+=y1;
        vars->gases[i].qz+=z1;
        //  set velocities
        vars->gases[i+Nof_around_gas].px=vars->gases[i].px+vx2;
        vars->gases[i+Nof_around_gas].py=vars->gases[i].py+vy2;
        vars->gases[i+Nof_around_gas].pz=vars->gases[i].pz+vz2;
        vars->gases[i].px+=vx1;
        vars->gases[i].py+=vy1;
        vars->gases[i].pz+=vz1;
        
        vars->gases[i].mass=14.01;
        vars->gases[i+Nof_around_gas].mass=14.01;
    }
}

void 
MD::makeDiatomicProp_out(int i){
    if(vars->gases[i+Nof_around_gas].mass!=0){
                //  averaged position
        vars->gases[i].qx=(vars->gases[i].qx+vars->gases[i+Nof_around_gas].qx)*0.5;
        vars->gases[i].qy=(vars->gases[i].qy+vars->gases[i+Nof_around_gas].qy)*0.5;
        vars->gases[i].qz=(vars->gases[i].qz+vars->gases[i+Nof_around_gas].qz)*0.5;
        //  averaged velocity
        vars->gases[i].px=(vars->gases[i].px+vars->gases[i+Nof_around_gas].px)*0.5;
        vars->gases[i].py=(vars->gases[i].py+vars->gases[i+Nof_around_gas].py)*0.5;
        vars->gases[i].pz=(vars->gases[i].pz+vars->gases[i+Nof_around_gas].pz)*0.5;
        
        vars->gases[i].mass=28.02;
        vars->gases[i+Nof_around_gas].mass=0;
        
        vars->gases[i+Nof_around_gas].px=0;
        vars->gases[i+Nof_around_gas].py=0;
        vars->gases[i+Nof_around_gas].pz=0;
    }
}