#include"macro.h"
#include <stdlib.h>
#include<string.h>
#include <stdio.h>
#include"initialize.h"
#include<gsl/gsl_rng.h>
#include <math.h>
#include "birth_pulsars.h"
#include "detection.h"
#include "evolution.h"
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf_gamma.h>
#include"galac_pot.h"

void distrib_vinit(void *params){

    struct func_params *part= (struct func_params*)params;
    long np;
    double vx,vy,vz;
    double age_pulsar_yr;
    const double yr_sec=365*24*3600;

    for(np=0;np<part->Npulsars;np++){
       
       age_pulsar_yr    =  part->age_pulsar[np]/yr_sec;
       if (age_pulsar_yr<3*1e6) part->sigma_v=part->v_young;
       else part->sigma_v=part->v_old;

       vx = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v); //km/s
       vy = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
       vz = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);

       part->vx0[np]=vx;part->vx[np]=vx;
       part->vy0[np]=vy;part->vy[np]=vy;
       part->vz0[np]=vz;part->vz[np]=vz;

    }

}

double phi_tot(void *params,long np){

      struct func_params *part= (struct func_params*)params;
      const double kpc2km=3.0856775807e16;
      double x=part->x[np]*kpc2km;double y=part->y[np]*kpc2km;double z=part->z[np]*kpc2km;
      double R=sqrt(sq(x)+sq(y));double r=sqrt(sq(x)+sq(y)+sq(z));
      double a1=0,b1=0.277*kpc2km,M1=1.12e10*MSUN;
      double a2=3.7*kpc2km,b2=0.2*kpc2km,M2=8.07e10*MSUN;
      double rc=6*kpc2km,Mc=5e10*MSUN;
      double phi;double phi_1;double phi_2;double phi_h;
      phi_1=(G_grav*M1)/(sqrt(sq(R)+sq(a1+sqrt(sq(z)+sq(b1)))));
      phi_2=(G_grav*M2)/(sqrt(sq(R)+sq(a2+sqrt(sq(z)+sq(b2)))));
      phi_h=-(G_grav*Mc/(rc))*(0.5*log(1+(sq(r)/sq(rc)))+(rc/r)*atan(r/rc));
      phi=phi_1+phi_2+phi_h;
      return phi;

}

double tot_energy(void *params,long np){

      struct func_params *part= (struct func_params*)params;
      double vx,vy,vz;
      vx=part->vx[np];vy=part->vy[np];vz=part->vz[np];
      double phi=phi_tot(part,np);
      double E_tot=0.5*(sq(vx)+sq(vy)+sq(vz))+phi;
      return E_tot;

}

void phi_1(void *params,double *k,double step,double phi1[3],long np){

       struct func_params *part= (struct func_params*)params;
       const double kpc2km=3.0856775807e16;
       double a1=0.0,b1=0.277*kpc2km,M1=1.12e10*MSUN;
       double x=part->x[np]*kpc2km;double y=part->y[np]*kpc2km;double z=part->z[np]*kpc2km;
       double x_rk=x+step*k[0];
       double y_rk=y+step*k[1];
       double z_rk=z+step*k[2];
       double R=sqrt(sq(x)+sq(y));//,Rx=sqrt(sq(x_rk)+sq(y)),Ry=sqrt(sq(x)+sq(y_rk));
       //double phi1_x=(G_grav*M1*x_rk*Rx)/(sqrt((sq(x_rk)+sq(y)+sq(a1+sqrt(sq(z)+b1))))*sqrt(sq(x_rk)+sq(y)+sq(z)));
       double phi1_x=-(2.0*G_grav*M1*x_rk)/(pow(sq(x_rk)+sq(y)+sq(a1+sqrt(sq(z)+sq(b1))),1.5));
       double phi1_y=-(2.0*G_grav*M1*y_rk)/(pow(sq(x)+sq(y_rk)+sq(a1+sqrt(sq(z)+sq(b1))),1.5));
       double phi1_z=-(G_grav*M1*z_rk)/(2.0*(pow(sq(R)+sq(a1+sqrt(sq(z_rk)+sq(b1))),1.5)*sqrt(sq(z_rk)+sq(b1))));
       //double phi1_y=(G_grav*M1*y_rk*Ry)/(sqrt((sq(x)+sq(y_rk)+sq(a1+sqrt(sq(z)+b1))))*sqrt(sq(x)+sq(y_rk)+sq(z)));
       //double phi1_z=(G_grav*M1*z_rk*R)/(sqrt((sq(x)+sq(y)+sq(a1+sqrt(sq(z_rk)+b1))))*sqrt(sq(x)+sq(y)+sq(z_rk)));
       phi1[0]=phi1_x;phi1[1]=phi1_y;phi1[2]=phi1_z;
}

void phi_2(void *params,double *k,double step,double phi2[3],long np){

       struct func_params *part= (struct func_params*)params;
       const double kpc2km=3.0856775807e16;
       double a2=3.7*kpc2km,b2=0.2*kpc2km,M2=8.07e10*MSUN;
       double x=part->x[np]*kpc2km;double y=part->y[np]*kpc2km;double z=part->z[np]*kpc2km;
       double x_rk=x+step*k[0];
       double y_rk=y+step*k[1];
       double z_rk=z+step*k[2];
       double R=sqrt(sq(x)+sq(y));//,Rx=sqrt(sq(x_rk)+sq(y)),Ry=sqrt(sq(x)+sq(y_rk));
       //double phi2_x=(G_grav*M2*x_rk*Rx)/(sqrt((sq(x_rk)+sq(y)+sq(a2+sqrt(sq(z)+b2))))*sqrt(sq(x_rk)+sq(y)+sq(z)));
       //double phi2_y=(G_grav*M2*y_rk*Ry)/(sqrt((sq(x)+sq(y_rk)+sq(a2+sqrt(sq(z)+b2))))*sqrt(sq(x)+sq(y_rk)+sq(z)));
       double phi2_x=-(2.0*G_grav*M2*x_rk)/(pow(sq(x_rk)+sq(y)+sq(a2+sqrt(sq(z)+sq(b2))),1.5));
       double phi2_y=-(2.0*G_grav*M2*y_rk)/(pow(sq(x)+sq(y_rk)+sq(a2+sqrt(sq(z)+sq(b2))),1.5));
       double phi2_z=-(G_grav*M2*z_rk)/(2.0*(pow(sq(R)+sq(a2+sqrt(sq(z_rk)+sq(b2))),1.5)*sqrt(sq(z_rk)+sq(b2))));
       //double phi2_z=(G_grav*M2*z_rk*R)/(sqrt((sq(x)+sq(y)+sq(a2+sqrt(sq(z_rk)+b2))))*sqrt(sq(x)+sq(y)+sq(z_rk)));
       phi2[0]=phi2_x;phi2[1]=phi2_y;phi2[2]=phi2_z;
}

void phih(void *params,double *k,double step,double phih[3],long np){

       struct func_params *part= (struct func_params*)params;
       const double kpc2km=3.0856775807e16;
       double rc=6*kpc2km,Mc=5e10*MSUN;
       double x=part->x[np]*kpc2km;double y=part->y[np]*kpc2km;double z=part->z[np]*kpc2km;
       double x_rk=x+step*k[0];
       double y_rk=y+step*k[1];
       double z_rk=z+step*k[2];
       double rx=sqrt(sq(x+step*k[0])+sq(y)+sq(z));
       double ry=sqrt(sq(x)+sq(y+step*k[1])+sq(z));
       double rz=sqrt(sq(x)+sq(y)+sq(z+step*k[2]));
       double phih_x=-(G_grav*Mc/rc)*((x_rk/(sq(rc)+sq(rx)))-((rc/(sq(rx)*rx))*atan(rx/rc)*x_rk)+(x_rk/sq(rx))*(1/(1+sq(rx/rc))));
       double phih_y=-(G_grav*Mc/rc)*((y_rk/(sq(rc)+sq(ry)))-((rc/(sq(ry)*ry))*atan(ry/rc)*y_rk)+(y_rk/sq(ry))*(1/(1+sq(ry/rc))));
       double phih_z=-(G_grav*Mc/rc)*((z_rk/(sq(rc)+sq(rz)))-((rc/(sq(rz)*rz))*atan(rz/rc)*z_rk)+(z_rk/sq(rz))*(1/(1+sq(rz/rc))));
       //double phih_x=(G_grav*Mc/(2.0*rc))*((2.0*rx/(sq(rc)+sq(rx)))-(rc/sq(rx))*atan(rx/rc)+((rc/rx)*(1.0/(1.0+sq(rx/rc)))))*((*part->x+step*k[0])/rx);
       //double phih_y=(G_grav*Mc/(2.0*rc))*((2.0*ry/(sq(rc)+sq(ry)))-(rc/sq(ry))*atan(ry/rc)+((rc/ry)*(1.0/(1.0+sq(ry/rc)))))*((*part->y+step*k[1])/ry);
       //double phih_z=(G_grav*Mc/(2.0*rc))*((2.0*rz/(sq(rc)+sq(rz)))-(rc/sq(rz))*atan(rz/rc)+((rc/rz)*(1.0/(1.0+sq(rz/rc)))))*((*part->z+step*k[2])/rz);
       phih[0]=phih_x;phih[1]=phih_y;phih[2]=phih_z;

}

void evol_galac_pot_RK4(void *params){

    struct func_params *part= (struct func_params*)params;
    const double yr_sec=365*24*3600;long np;
    double step=10e4*yr_sec;
    FILE *file=NULL;
    const double kpc2km=3.0856775807e16;
    file=fopen("x_y_err_RK4.txt","w+");
    for(np=0;np<part->Npulsars;np++){

       double E0=tot_energy(part,np);
       for(double t=TMILKY*yr_sec-part->age_pulsar[np];t<TMILKY*yr_sec;t+=step){

	 double vx,vy,vz;
	 double x,y,z;
	 double gphih_1[3],gphi_1_1[3],gphi_2_1[3];
	 double gphih_2[3],gphi_1_2[3],gphi_2_2[3];
	 double gphih_3[3],gphi_1_3[3],gphi_2_3[3];
	 double gphih_4[3],gphi_1_4[3],gphi_2_4[3];
         double k1[3];double k2[3];double k3[3];double k4[3];double k0[3];
	 k0[0]=0;k0[1]=0;k0[2]=0;
         int i;
	 double kx_1,kx_2,kx_3,kx_4,ky_1,ky_2,ky_3,ky_4,kz_1,kz_2,kz_3,kz_4; 
	 phih(part,k0,0,gphih_1,np);phi_1(part,k0,0,gphi_1_1,np);phi_2(part,k0,0,gphi_2_1,np);
	 for(i=0;i<3;i++){

		 k1[i]=-gphih_1[i]-gphi_1_1[i]-gphi_2_1[i];

	 }
	 phih(part,k1,step*0.5,gphih_2,np);phi_1(part,k1,step*0.5,gphi_1_2,np);phi_2(part,k1,step*0.5,gphi_2_2,np);
	 for(i=0;i<3;i++){

                 k2[i]=-gphih_2[i]-gphi_1_2[i]-gphi_2_2[i];

         }
	 phih(part,k2,0.5*step,gphih_3,np);phi_1(part,k2,0.5*step,gphi_1_3,np);phi_2(part,k2,0.5*step,gphi_2_3,np);
	 for(i=0;i<3;i++){

                 k3[i]=-gphih_3[i]-gphi_1_3[i]-gphi_2_3[i];

         }
	 phih(part,k3,step,gphih_4,np);phi_1(part,k3,step,gphi_1_4,np);phi_2(part,k3,step,gphi_2_4,np);
	 for(i=0;i<3;i++){

                 k4[i]=-gphih_4[i]-gphi_1_4[i]-gphi_2_4[i];

         }

	 vx=part->vx[np]+(step/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
	 vy=part->vy[np]+(step/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
	 vz=part->vz[np]+(step/6.0)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);

         kx_1=part->vx[np];ky_1=part->vy[np];kz_1=part->vz[np];
	 kx_2=part->vx[np]+step*0.5*k1[0];ky_2=part->vy[np]+step*0.5*k1[1];kz_2=part->vz[np]+step*0.5*k1[2];
	 kx_3=part->vx[np]+step*0.5*k2[0];ky_3=part->vy[np]+step*0.5*k2[1];kz_3=part->vz[np]+step*0.5*k2[2];
	 kx_4=part->vx[np]+step*k3[0];ky_4=part->vy[np]+step*k3[1];kz_4=part->vz[np]+step*k3[2];

	 x=part->x[np]+(step/6.0)*((kx_1+2*kx_2+2*kx_3+kx_4)/kpc2km);
	 y=part->y[np]+(step/6.0)*((ky_1+2*ky_2+2*ky_3+ky_4)/kpc2km);
	 z=part->z[np]+(step/6.0)*((kz_1+2*kz_2+2*kz_3+kz_4)/kpc2km);

	 part->vx[np]=vx;part->vy[np]=vy;part->vz[np]=vz;
	 part->x[np]=x;part->y[np]=y;part->z[np]=z;

       }

       part->x_s[np]= 8.5-part->x[np]; // shift center of the Galaxy to the Sun
       part->y_s[np]= part->y[np];
       part->z_s[np]= 0.015-part->z[np];

       double Ef=tot_energy(part,np);
       part->err_rel_g[np]=fabs((E0-Ef)/E0)*100;
       fprintf(file,"%e %e %e \n",part->x[np],part->y[np],part->err_rel_g[np]);
}

    fclose(file);

}

void evol_galac_pot_verlet(void *params){

    struct func_params *part= (struct func_params*)params;
    const double yr_sec=365*24*3600;long np;
    double step=10e2*yr_sec;
    FILE *file=NULL;
    const double kpc2km=3.0856775807e16;
    file=fopen("x_y_err_verlet.txt","w+");
    for(np=0;np<part->Npulsars;np++){

       double E0=tot_energy(part,np);
       for(double t=TMILKY*yr_sec-part->age_pulsar[np];t<TMILKY*yr_sec;t+=step){

         double gphih[3],gphi_1[3],gphi_2[3];
         double k0[3];double grad_phi[3];
         k0[0]=0;k0[1]=0;k0[2]=0;
	 int i;

	 //First shift of the space coordinates
	 part->x[np]=part->x[np]+(part->vx[np]*step*0.5)/kpc2km;
	 part->y[np]=part->y[np]+(part->vy[np]*step*0.5)/kpc2km;
	 part->z[np]=part->z[np]+(part->vz[np]*step*0.5)/kpc2km;

	 //Computation of the gradient with the new space coordinates
	 phih(part,k0,0,gphih,np);phi_1(part,k0,0,gphi_1,np);phi_2(part,k0,0,gphi_2,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphih[i]+gphi_1[i]+gphi_2[i];

         }

	 //Shift of the velocities 
	 part->vx[np]=part->vx[np]-grad_phi[0]*step;
	 part->vy[np]=part->vy[np]-grad_phi[1]*step;
	 part->vz[np]=part->vz[np]-grad_phi[2]*step;

	 //Second shift of the coordinates 
	 part->x[np]=part->x[np]+(part->vx[np]*step*0.5)/kpc2km;
         part->y[np]=part->y[np]+(part->vy[np]*step*0.5)/kpc2km;
         part->z[np]=part->z[np]+(part->vz[np]*step*0.5)/kpc2km;


       }

     part->x_s[np]= 8.5-part->x[np]; // shift center of the Galaxy to the Sun
     part->y_s[np]= part->y[np];
     part->z_s[np]= 0.015-part->z[np];

     double Ef=tot_energy(part,np);
     part->err_rel_g[np]=fabs((E0-Ef)/E0)*100;
     fprintf(file,"%e %e %e \n",part->x[np],part->y[np],part->err_rel_g[np]);

    }
    
  fclose(file);

}

void evol_galac_pot_yoshi(void *params){

    struct func_params *part= (struct func_params*)params;
    const double yr_sec=365*24*3600;long np;
    double step=10e2*yr_sec;
    FILE *file=NULL;
    const double kpc2km=3.0856775807e16;
    double w0=-(pow(2,1.5))/(2-pow(2,1.5));double w1=(1/(2-pow(2,1.5)));
    double c1=w1/2.0;double c4=c1;double c2=(w0+w1)/2.0;double c3=c2;
    double d1=w1;double d3=w1;double d2=w0;
    file=fopen("x_y_err_verlet.txt","w+");
    for(np=0;np<part->Npulsars;np++){

       double E0=tot_energy(part,np);
       for(double t=TMILKY*yr_sec-part->age_pulsar[np];t<TMILKY*yr_sec;t+=step){

         double gphih[3],gphi_1[3],gphi_2[3];
	 double gphih_2[3],gphi_1_2[3],gphi_2_2[3];
	 double gphih_3[3],gphi_1_3[3],gphi_2_3[3];
         double k0[3];double grad_phi[3];
         k0[0]=0;k0[1]=0;k0[2]=0;
         int i;
	 // first shift coordinates
	 part->x[np]=part->x[np]+(c1*part->vx[np]*step)/kpc2km;
	 part->y[np]=part->y[np]+(c1*part->vy[np]*step)/kpc2km;
	 part->z[np]=part->z[np]+(c1*part->vz[np]*step)/kpc2km;
	 
	 // first shift velocities 
	 phih(part,k0,0,gphih,np);phi_1(part,k0,0,gphi_1,np);phi_2(part,k0,0,gphi_2,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphih[i]+gphi_1[i]+gphi_2[i];

         }
	 part->vx[np]=part->vx[np]+d1*grad_phi[0]*step;
	 part->vy[np]=part->vy[np]+d1*grad_phi[1]*step;
	 part->vz[np]=part->vz[np]+d1*grad_phi[2]*step;

	 //second shift coordinates
	 part->x[np]=part->x[np]+(c2*part->vx[np]*step)/kpc2km;
         part->y[np]=part->y[np]+(c2*part->vy[np]*step)/kpc2km;
         part->z[np]=part->z[np]+(c2*part->vz[np]*step)/kpc2km;

	 //second shift velocities
	 phih(part,k0,0,gphih_2,np);phi_1(part,k0,0,gphi_1_2,np);phi_2(part,k0,0,gphi_2_2,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphih_2[i]+gphi_1_2[i]+gphi_2_2[i];

         }
         part->vx[np]=part->vx[np]+d2*grad_phi[0]*step;
         part->vy[np]=part->vy[np]+d2*grad_phi[1]*step;
         part->vz[np]=part->vz[np]+d2*grad_phi[2]*step;

	 //third shift coordinates
	 part->x[np]=part->x[np]+(c3*part->vx[np]*step)/kpc2km;
         part->y[np]=part->y[np]+(c3*part->vy[np]*step)/kpc2km;
         part->z[np]=part->z[np]+(c3*part->vz[np]*step)/kpc2km;

	 //third shift of velocities
	 phih(part,k0,0,gphih_3,np);phi_1(part,k0,0,gphi_1_3,np);phi_2(part,k0,0,gphi_2_3,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphih_3[i]+gphi_1_3[i]+gphi_2_3[i];

         }
         part->vx[np]=part->vx[np]+d3*grad_phi[0]*step;
         part->vy[np]=part->vy[np]+d3*grad_phi[1]*step;
         part->vz[np]=part->vz[np]+d3*grad_phi[2]*step;

	 //fourt shift of coordinates
	 part->x[np]=part->x[np]+(c4*part->vx[np]*step)/kpc2km;
         part->y[np]=part->y[np]+(c4*part->vy[np]*step)/kpc2km;
         part->z[np]=part->z[np]+(c4*part->vz[np]*step)/kpc2km;

	 //fourth shift of velocities : nothing to do
	 


       }

       part->x_s[np]= 8.5-part->x[np]; // shift center of the Galaxy to the Sun
       part->y_s[np]= part->y[np];
       part->z_s[np]= 0.015-part->z[np];

       double Ef=tot_energy(part,np);
       part->err_rel_g[np]=fabs((E0-Ef)/E0)*100;
       fprintf(file,"%e %e %e \n",part->x[np],part->y[np],part->err_rel_g[np]);

       }

    fclose(file);
}

void send_data(void *params){

     struct func_params *part= (struct func_params*)params;
     long np;
     FILE *file=NULL;
     file=fopen("data_for_py.txt","w+");
     for(np=0;np<part->Npulsars;np++){

	 double R=sqrt(sq(part->x[np])+sq(part->y[np]));
	 double phi=atan(part->y[np]/part->x[np]);
	 double vr=(part->x[np]*part->vx[np]+part->y[np]*part->vy[np])/(sqrt(sq(part->x[np])+sq(part->y[np])));
	 double vphi=-(1.0/(1+sq(part->y[np]/part->x[np])))*((part->y[np]*(part->vx[np]/kpc2km)/(2.0*sq(part->x[np])))+(part->x[np]*(part->vy[np]/kpc2km)/(2.0*sq(part->y[np]))));
	 fprintf(file,"%e|%e|%e|%e|%e|%e|\n",R,vr,vphi,part->z[np],part->vz[np],phi);
     }
     fclose(file);


}
