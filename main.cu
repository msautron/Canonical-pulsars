#define GNU_SOURCE
//#define _GNU_SOURCE
#include <stdlib.h>
#include"initialize.h"
#include"macro.h"
#include"birth_pulsars.h"
#include"evolution.h"
#include"galac_pot.h"
#include"detection.h"
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include "cn.h"
#include <curand_kernel.h>
#include <math_constants.h>
#define TAILLE_MAX 1000 // Tableau de taille 1000

#define NTHREAD (250)
#define NBLOCK (500)
#define NPOINTS (NBLOCK*NTHREAD)

__device__ double phi_tot(double *gx,double *gy,double *gz,int np){ //Compute the gravitational potential felt by a pulsar

      const double kpc2km=3.0856775807e16;
      double L0=1.0;
      double x=gx[np];double y=gy[np];double z=gz[np];
      double R=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km));double r=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(z*L0*kpc2km));
      double a1=0,b1=0.267*kpc2km,M1=1.02e10*MSUN;
      double a2=4.4*kpc2km,b2=0.308*kpc2km,M2=6.50535e10*MSUN;
      //double rc=6*kpc2km,Mc=5e10*MSUN;
      double Mh=12474*2.325*1e7*MSUN;double ah=7.7*kpc2km;
      double phi;
      double phi_1;
      double phi_2;
      double phi_k;
      double phi_NFW;
      phi_1=-(G_grav*M1)/(sqrt(sq(R)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1)))));
      phi_2=-(G_grav*M2)/(sqrt(sq(R)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2)))));
      phi_NFW=-(G_grav*Mh/(r))*log(1.0+(r/ah));
      phi_k=-(G_grav*4e6*MSUN)/r;
      phi=phi_1+phi_2+phi_NFW+phi_k;
      return phi;

}

__device__ double tot_energy(double *gx,double *gy,double *gz,double *gvx,double *gvy,double *gvz,int np){ //Compute the sum of Kinetic and potential energy of a pulsar

      double vx,vy,vz;
      double v0=100.0;
      vx=gvx[np]*v0;vy=gvy[np]*v0;vz=gvz[np]*v0;
      double phi=phi_tot(gx,gy,gz,np);
      double E_tot=0.5*(sq(vx)+sq(vy)+sq(vz))+phi;
      return E_tot;

}

__device__ void phi_1(double *gx,double *gy,double *gz,double phi1[3],int np){ //Potential for the bulge of the Milky Way

       const double kpc2km=3.0856775807e16;
       double a1=0.0,b1=0.267*kpc2km,M1=1.02e10*MSUN;
       double L0=1.0;double v0=100.0;
       double T0=(L0*kpc2km)/v0;
       double x=gx[np];double y=gy[np];double z=gz[np];
       double R=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km));
       double phi1_x=-((G_grav*M1*x)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1))),1.5)))*(sq(T0)/L0);
       double phi1_y=-((G_grav*M1*y)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1))),1.5)))*(sq(T0)/L0);
       double phi1_z=-((G_grav*M1*z)/((pow(sq(R)+sq(a1+sqrt(sq(z*L0*kpc2km)+sq(b1))),1.5)*sqrt(sq(z*L0*kpc2km)+sq(b1)))))*(sq(T0)/L0)*(a1+sqrt(sq(z*L0*kpc2km)+sq(b1)));
       phi1[0]=phi1_x;phi1[1]=phi1_y;phi1[2]=phi1_z;
}

__device__ void phi_2(double *gx,double *gy,double *gz,double phi2[3],int np){ //Potential for the disk of the Milky Way

       const double kpc2km=3.0856775807e16;
       double a2=4.4*kpc2km,b2=0.308*kpc2km,M2=6.50535e10*MSUN;
       double L0=1.0;double v0=100.0;
       double T0=(L0*kpc2km)/v0;
       double x=gx[np];double y=gy[np];double z=gz[np];
       double R=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km));
       double phi2_x=-((G_grav*M2*x)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2))),1.5)))*(sq(T0)/L0);
       double phi2_y=-((G_grav*M2*y)/(pow(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2))),1.5)))*(sq(T0)/L0);
       double phi2_z=-((G_grav*M2*z)/((pow(sq(R)+sq(a2+sqrt(sq(z*L0*kpc2km)+sq(b2))),1.5)*sqrt(sq(z*L0*kpc2km)+sq(b2)))))*(sq(T0)/L0)*(a2+sqrt(sq(z*L0*kpc2km)+sq(b2)));
       phi2[0]=phi2_x;phi2[1]=phi2_y;phi2[2]=phi2_z;
}

__device__ void phi_NFW(double *gx,double *gy,double *gz,double phi_NFW[3],int np){ //Potential for the dark matter halo of the galaxy

        const double kpc2km=3.0856775807e16;
        double Mh=12474*2.325*1e7*MSUN;double ah=7.7*kpc2km;double L0=1.0;double v0=100.0;
        double T0=(L0*kpc2km)/v0;
        double x=gx[np];double y=gy[np];double z=gz[np];
        double r=sqrt(sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(z*L0*kpc2km));
        double phi_x=-(G_grav*Mh*x/(sq(r)))*(((1.0/r)*log(1.0+(r/ah)))-(1.0/(ah+r)))*(sq(T0)/L0);
        double phi_y=-(G_grav*Mh*y/(sq(r)))*(((1.0/r)*log(1.0+(r/ah)))-(1.0/(ah+r)))*(sq(T0)/L0);
        double phi_z=-(G_grav*Mh*z/(sq(r)))*(((1.0/r)*log(1.0+(r/ah)))-(1.0/(ah+r)))*(sq(T0)/L0);
        phi_NFW[0]=phi_x;phi_NFW[1]=phi_y;phi_NFW[2]=phi_z;

}

__device__ void phi_kepler(double *gx,double *gy,double *gz,double phi_k[3],int np){ //Potential for the supermassive black hole at the center of the Milky Way

        const double kpc2km=3.0856775807e16;
        double Mc=4e6*MSUN;double L0=1.0;double v0=100.0;
        double T0=(L0*kpc2km)/v0;
        double x=gx[np];double y=gy[np];double z=gz[np];
        double r_2=sq(x*L0*kpc2km)+sq(y*L0*kpc2km)+sq(z*L0*kpc2km);
        double phi_x=-(G_grav*Mc*x*sq(T0))/(pow(r_2,1.5)*L0);
        double phi_y=-(G_grav*Mc*y*sq(T0))/(pow(r_2,1.5)*L0);
        double phi_z=-(G_grav*Mc*z*sq(T0))/(pow(r_2,1.5)*L0);
        phi_k[0]=phi_x;phi_k[1]=phi_y;phi_k[2]=phi_z;

}

__global__ void evol_galac_PEFRL(double *gx,double *gy,double *gz, double *gvx,double *gvy,double *gvz,double *gage_pulsar,int *gsize,double *error_rel,double *gtmilky,double *ggl,double *ggb,double *gDM,double *gperiod,double *gw_r_fast,double *gw_r_pmps,double *galpha,double *g_nomegax,double *g_nomegay,double *g_nomegaz,double *gxi,double *grho,long seed,long npulsar,double *gnx,double *gny,double *gnz){ //PEFRL integration scheme + computation of DM and geometry of emission

    const double yr_sec=365*24*3600;int np;
    double step;
    const double kpc2km=3.0856775807e16;
    double L0=1.0;double v0=100.0;
    double T0=(L0*kpc2km)/v0;
    double step_wd;
    double xsi=0.1786178958448091;double lambda=-0.2123418310626054;
    double chi=-0.6626458266981849e-01;
    int idx=threadIdx.x+blockIdx.x*NTHREAD;
    long total_threads = gridDim.x*NTHREAD;
    long size=(npulsar + total_threads - 1) / total_threads;
    //int size=*gsize;
    long start=idx*size;//long end=start+size;
    long end=min(npulsar,start+size);
    double x_s;double y_s;double dist;double z_s;
    //printf("Bonjour I'm thread %d in block %d, with absolute %d\n",threadIdx.x,blockIdx.x,idx);
    for(np=start;np<end;np++){

       //printf("Computation for pulsar nb=%d in thread %d in block %d, with absolute %d\n",np,threadIdx.x,blockIdx.x,idx);
       //No dimension for the variables
       gx[np]=gx[np]/L0;gy[np]=gy[np]/L0;gz[np]=gz[np]/L0;
       gvx[np]=gvx[np]/v0;gvy[np]=gvy[np]/v0;gvz[np]=gvz[np]/v0;

       //Beginning of computation
       double E0=tot_energy(gx,gy,gz,gvx,gvy,gvz,np);
       for(double t=(*gtmilky)-gage_pulsar[np];t<(*gtmilky);t+=step){

	 if(gage_pulsar[np]<=1e2*yr_sec) {step=10*yr_sec;}
	 else if(gage_pulsar[np]>1e2*yr_sec && gage_pulsar[np]<=1e3*yr_sec) {step=25*yr_sec;}
	 else if(gage_pulsar[np]>1e3*yr_sec && gage_pulsar[np]<=1e4*yr_sec) {step=5e1*yr_sec;}
	 else if(gage_pulsar[np]>1e4*yr_sec && gage_pulsar[np]<=1e5*yr_sec) {step=5e2*yr_sec;}
	 else if(gage_pulsar[np]>1e5*yr_sec && gage_pulsar[np]<=1e6*yr_sec) {step=2.5e3*yr_sec;}
	 else if(gage_pulsar[np]>1e6*yr_sec && gage_pulsar[np]<=1e7*yr_sec) {step=2.5e4*yr_sec;}
	 else if(gage_pulsar[np]>1e7*yr_sec && gage_pulsar[np]<=1e8*yr_sec) {step=1e5*yr_sec;}
	 else if(gage_pulsar[np]>1e8*yr_sec) {step=1e5*yr_sec;}

	 step_wd=step/T0;

         double gphi_1[3];
         double gphi_2[3];
         double gphi_NFW[3];
         double grad_phi[3];
	 double gphi_k[3];
         int i;

         //First shift of space coordinates
         gx[np]+=gvx[np]*step_wd*xsi;
         gy[np]+=gvy[np]*step_wd*xsi;
         gz[np]+=gvz[np]*step_wd*xsi;

         //Computation of the gradient with the new space coordinates
         phi_1(gx,gy,gz,gphi_1,np);
         phi_2(gx,gy,gz,gphi_2,np);
         phi_NFW(gx,gy,gz,gphi_NFW,np);
	 phi_kepler(gx,gy,gz,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //Shift of the velocities
         gvx[np]+=grad_phi[0]*step_wd*0.5*(1-2*lambda);
         gvy[np]+=grad_phi[1]*step_wd*0.5*(1-2*lambda);
         gvz[np]+=grad_phi[2]*step_wd*0.5*(1-2*lambda);

         //second shift of space coordinates
         gx[np]+=gvx[np]*step_wd*chi;
         gy[np]+=gvy[np]*step_wd*chi;
         gz[np]+=gvz[np]*step_wd*chi;

         //Computation of the gradient with the new space coordinates
	 phi_1(gx,gy,gz,gphi_1,np);
         phi_2(gx,gy,gz,gphi_2,np);
         phi_NFW(gx,gy,gz,gphi_NFW,np);
         phi_kepler(gx,gy,gz,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //Second shift of velocities
         gvx[np]+=grad_phi[0]*step_wd*lambda;
         gvy[np]+=grad_phi[1]*step_wd*lambda;
         gvz[np]+=grad_phi[2]*step_wd*lambda;

         //third shift of space coordinates
         gx[np]+=gvx[np]*step_wd*(1-2*(chi+xsi));
         gy[np]+=gvy[np]*step_wd*(1-2*(chi+xsi));
         gz[np]+=gvz[np]*step_wd*(1-2*(chi+xsi));

         //Computation of the gradient with the new space coordinates
	 phi_1(gx,gy,gz,gphi_1,np);
         phi_2(gx,gy,gz,gphi_2,np);
         phi_NFW(gx,gy,gz,gphi_NFW,np);
         phi_kepler(gx,gy,gz,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //third shift of velocities
         gvx[np]+=grad_phi[0]*step_wd*lambda;
         gvy[np]+=grad_phi[1]*step_wd*lambda;
         gvz[np]+=grad_phi[2]*step_wd*lambda;

         //fourth shift of space coordinates
         gx[np]+=gvx[np]*step_wd*chi;
         gy[np]+=gvy[np]*step_wd*chi;
         gz[np]+=gvz[np]*step_wd*chi;

         //Computation of the gradient with the new space coordinates
	 phi_1(gx,gy,gz,gphi_1,np);
         phi_2(gx,gy,gz,gphi_2,np);
         phi_NFW(gx,gy,gz,gphi_NFW,np);
         phi_kepler(gx,gy,gz,gphi_k,np);
         for(i=0;i<3;i++){

                 grad_phi[i]=gphi_NFW[i]+gphi_1[i]+gphi_2[i]+gphi_k[i];

            }

         //fourth shift of velocities
         gvx[np]+=grad_phi[0]*step_wd*0.5*(1-2*lambda);
         gvy[np]+=grad_phi[1]*step_wd*0.5*(1-2*lambda);
         gvz[np]+=grad_phi[2]*step_wd*0.5*(1-2*lambda);

         //final shift space coord
         gx[np]+=gvx[np]*step_wd*xsi;
         gy[np]+=gvy[np]*step_wd*xsi;
         gz[np]+=gvz[np]*step_wd*xsi;
       }

      //Save values in a file + computation energy and error
      double Ef=tot_energy(gx,gy,gz,gvx,gvy,gvz,np);
      error_rel[np]=fabs((E0-Ef)/E0)*100;
      gx[np]=gx[np]*L0;gy[np]=gy[np]*L0;gz[np]=gz[np]*L0;
      gvx[np]=gvx[np]*v0;gvy[np]=gvy[np]*v0;gvz[np]=gvz[np]*v0;

      x_s= gx[np]; // shift center of the Galaxy to the Sun
      y_s= gy[np]-8.5;
      z_s= gz[np]-0.015;
      dist= sqrt(sq(x_s)+sq(y_s)+sq(z_s));
      double r= sqrt(sq(x_s)+sq(y_s));
      double gbr;double glr;double RAD_2=180/M_PI;
      if(dist<1e-15){
      		gbr=0;
                }
      else{
                gbr=asin(z_s/dist);
                }
      ggb[np]=gbr*RAD_2;

      if(r<1e-15){
                glr=0;
                  }
      else{
      		if(x_s>=0){
                      glr=acos(-y_s/r);
                      }
                else{
                      glr=acos(y_s/r)+M_PI;
                    }
          }

      ggl[np]=glr*RAD_2;
      //Compute the geometry of emission
      double h_em=2e5;
      double alpha;double rho;double norm;double n[3];double xi;double ratio;
      /*double cos_theta_g;
      double phi_g;
      curandState state;
      curand_init(seed, idx, 0, &state);
      cos_theta_g=2*curand_uniform(&state)-1;
      phi_g=2*M_PI*curand_uniform(&state);
      g_nomegax[np]=((1-pow(cos_theta_g,2.0))*cosf(phi_g));
      g_nomegay[np]=((1-pow(cos_theta_g,2.0))*sinf(phi_g));
      g_nomegaz[np]=cos_theta_g;*/
      if (galpha[np]<=M_PI/2) alpha=galpha[np]; // just to make it easier to read    
      else if (galpha[np]>M_PI/2) alpha=M_PI-galpha[np];
      rho            = 3*sqrt(M_PI*h_em)/sqrt(2*gperiod[np]*SI_C); // 3*sqrt((PI*hem)/2*P*c) Equation 2 of Johnston et al. (2020)
      grho[np]=rho;

      //printf("rho : %e\n",180*rho/M_PI);

      /* line of sight vector */
      norm=sqrt(sq(gx[np])+sq(gy[np]-8.5)+sq(gz[np]-0.015));
      n[0]=(gx[np])/norm;
      n[1]=(gy[np]-8.5)/norm;
      n[2]=(gz[np]-0.015)/norm;
      gnx[np]=n[0];gny[np]=n[1];gnz[np]=n[2];

      /* angle between vectors n and n_omega  */
      xi=acos(g_nomegax[np]*n[0]+g_nomegay[np]*n[1]+g_nomegaz[np]*n[2]);
      gxi[np]=xi;
      if (gxi[np]<=M_PI/2) xi=gxi[np]; // just to make it easier to read
      else if (gxi[np]>M_PI/2) xi=M_PI-gxi[np];



      /* calculates the width of the radio beam w_r, from eq 22 of our paper */
      if(fabs(alpha-xi)<= rho && alpha >= rho){
      		ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
                gw_r_fast[np]=2*acos(ratio);
      } else if(fabs(xi-(M_PI-alpha))<=rho && alpha >= rho){
                //alpha=M_PI-alpha;
                ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
                gw_r_fast[np]=2*acos(ratio);
                }

      //Compute the DM of each pulsar
      double DM_Host=100;int nt=1;int vbs=0;char dirname[10]="./";char text[10]="";int ndir=2;
      double tau_sc;
      double wint;
      double tau_samp_fast=0.05e-3; //tau_samp of FAST
      double e=1.6e-19;//electronic charge in C
      double m_e=9.1e-31; //electron mass in kg
      double f=1.374e9; //Frequency of survey fast/PMPS
      double deltaf_fast=0.24414e6; //bandwitdh of observation of fast
      double tau_samp_pmps=250e-6; //tau_samp of PMPS
      double delta_f_pmps=3000e3; //bandwidth of observation of PMPS
      /*double f=1400.0e6; //Frequency of survey SKA-mid
      double f_SKAlow=300e6; //Frequency of survey SKA-low
      double deltaf_fast=0.005e6; //bandwitdh of observation of SKA-low
      double tau_samp_pmps=6.4e-5; //tau_samp of SKA-mid
      double delta_f_pmps=0.009e6; //bandwidth of observation of SKA-mid
      double tau_samp_fast=0.1e-3; //tau_samp of SKA-low*/
      double tau_dm_fast,tau_dm_pmps;
      double DM_SI;
      wint=gw_r_fast[np]*gperiod[np]/(2*M_PI);
      dmdtau(ggl[np], ggb[np], dist*1e3, DM_Host, ndir, nt, vbs, dirname, text, &gDM[np],&tau_sc); //Compute the DM from the distance (Yao et al. (2017))
      DM_SI=gDM[np]*3.086e22;
      tau_dm_fast=(sq(e)*deltaf_fast*DM_SI)/(4*M_PI*SI_eps0*M_PI*m_e*SI_C*cube(f));
      tau_dm_pmps=(sq(e)*delta_f_pmps*DM_SI)/(4*M_PI*SI_eps0*M_PI*m_e*SI_C*cube(f));
      gw_r_fast[np]=sqrt(sq(wint)+sq(tau_samp_fast)+sq(tau_dm_fast)+sq(tau_sc))*((2*M_PI)/(gperiod[np]));
      gw_r_pmps[np]=sqrt(sq(wint)+sq(tau_samp_pmps)+sq(tau_dm_pmps)+sq(tau_sc))*((2*M_PI)/(gperiod[np]));
      //if (isnan(gw_r[np])) printf("w_int=%e, tau_dm=%e, tau_sc=%e, DM=%e, period=%e, gl=%e, gb=%e, dist=%e\n",wint,tau_dm,tau_sc,gDM[np],gperiod[np],ggl[np],ggb[np],dist);
    }
    printf("The thread %d in block %d, with absolute %d has finished\n",threadIdx.x,blockIdx.x,idx);
}

int main(int argc, char **argv){

	time_t start_time,end_time;
	time(&start_time);
	int nb_gpu=0;

	cudaError_t cudaStatus=cudaSetDevice(nb_gpu);
	if(cudaStatus!=cudaSuccess){
		fprintf(stderr, "Erreur lors de la sélection du GPU : %s\n", cudaGetErrorString(cudaStatus));
		return 1;
	}

	dim3 dimblock(NTHREAD);
        dim3 dimgrid(NBLOCK);

        struct func_params params;
    	params.T = gsl_rng_default;
    	params.r = gsl_rng_alloc (params.T);

	 time_t result = time(NULL);
         if(result != (time_t)(-1)) gsl_rng_set(params.r,(intmax_t)result);
             printf("# The current time is %s(%jd seconds since the Epoch)\n",
                  asctime(gmtime(&result)), (intmax_t)result);

        initialize(argc, argv, &params); //parameters initialization

	srand((unsigned)time(NULL));
	gsl_rng_set(params.r, time(NULL));

	printf("## birth_rate (yr) %ld \n",params.birth_rate2);
	printf("## sigma_P (s or log(s)) %e \n",params.sigma_p);
	printf("## sigma_B (logB) %e \n",params.sigma_b);
	printf("## Pmean (s) %e \n",params.p_mean);
	printf("## Bmean (T) %e \n",params.b_mean);
	printf("## alpha_d  %e \n",params.alpha_d);
	printf("## v_init (km.s-1) %ld \n",params.v_old);
	printf("\n");

         /* Calculates P, dotP and B */
	birth(&params); // generates the pulsar with an initial Pinit, Binit and age
        evolution(&params); //calculates P,Pdot and Edot
  
        // Distribution of the pulsars in the Galaxy
	//distrib_init(&params);    // initial distribution of the pulsars    
	distrib_init_2(&params); // initial distribution of the pulsars, when they are born in spiral arms
	distrib_vinit(&params); // Use when you want to run the simulation with the galactic potential
	
	double *gvx;double *gvy;double *gvz;double *gx;double *gy;double *gz;double *gage_pulsar;int *gsize;double *error_rel;double *gtmilky;double *ggl;double *ggb;double *gDM;double *gperiod;double *gw_r_fast;double *gw_r_pmps;double *galpha;double *g_nomegax;
	double *g_nomegay;double *g_nomegaz;double *gxi;double *grho;double *gnx;double *gny;double *gnz;
	double *x_bis;double *y_bis;double *z_bis;double *vx_bis;double *vy_bis;double *vz_bis;double *error_bis;double *gl_bis;double *gb_bis;double tmilky;double *DM_bis;double *w_r_fast_bis;double *w_r_pmps_bis;double *xi_bis;
	double *rho_bis;double *nx_bis;double *ny_bis;double *nz_bis;
	const double yr_sec=365*24*3600;tmilky=TMILKY*yr_sec;
	int size=params.Npulsars/(NBLOCK*NTHREAD);
        cudaMalloc((void**)&gvx,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gvy,sizeof(double)*params.Npulsars);
 	cudaMalloc((void**)&gvz,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gx,sizeof(double)*params.Npulsars);
 	cudaMalloc((void**)&gy,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gz,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&ggl,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&ggb,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gage_pulsar,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gsize,sizeof(int));
	cudaMalloc((void**)&error_rel,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gtmilky,sizeof(double));
	cudaMalloc((void**)&gDM,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gperiod,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gw_r_fast,sizeof(double)*params.Npulsars);
        cudaMalloc((void**)&gw_r_pmps,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&galpha,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&g_nomegax,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&g_nomegay,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&g_nomegaz,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gxi,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&grho,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gnx,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gny,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gnz,sizeof(double)*params.Npulsars);

	x_bis=(double *)calloc(params.Npulsars,sizeof(double));
	y_bis=(double *)calloc(params.Npulsars,sizeof(double));
	z_bis=(double *)calloc(params.Npulsars,sizeof(double));
	vx_bis=(double *)calloc(params.Npulsars,sizeof(double));
	vy_bis=(double *)calloc(params.Npulsars,sizeof(double));
	vz_bis=(double *)calloc(params.Npulsars,sizeof(double));
	error_bis=(double *)calloc(params.Npulsars,sizeof(double));
	gl_bis=(double *)calloc(params.Npulsars,sizeof(double));
	gb_bis=(double *)calloc(params.Npulsars,sizeof(double));
	DM_bis=(double *)calloc(params.Npulsars,sizeof(double));
	w_r_pmps_bis=(double *)calloc(params.Npulsars,sizeof(double));
        w_r_fast_bis=(double *)calloc(params.Npulsars,sizeof(double));
	xi_bis=(double *)calloc(params.Npulsars,sizeof(double));
	rho_bis=(double *)calloc(params.Npulsars,sizeof(double));
	nx_bis=(double *)calloc(params.Npulsars,sizeof(double));
	ny_bis=(double *)calloc(params.Npulsars,sizeof(double));
	nz_bis=(double *)calloc(params.Npulsars,sizeof(double));

        cudaMemcpy(gvx,params.vx,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gvy,params.vy,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gvz,params.vz,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gx,params.x,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gy,params.y,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gz,params.z,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(ggl,params.gl,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(ggb,params.gb,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gage_pulsar,params.age_pulsar,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gsize,&size,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(gtmilky,&tmilky,sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gDM,params.DM,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gperiod,params.period,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gw_r_pmps,params.w_r_pmps,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(gw_r_fast,params.w_r_fast,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(galpha,params.alpha,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(g_nomegax,params.n_omega_x,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(g_nomegay,params.n_omega_y,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(g_nomegaz,params.n_omega_z,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gxi,params.xi,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(grho,params.rho,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gnx,params.nx,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gny,params.ny,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gnz,params.nz,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);

        evol_galac_PEFRL<<<dimgrid,dimblock>>>(gx,gy,gz,gvx,gvy,gvz,gage_pulsar,gsize,error_rel,gtmilky,ggl,ggb,gDM,gperiod,gw_r_fast,gw_r_pmps,galpha,g_nomegax,g_nomegay,g_nomegaz,gxi,grho,time(NULL),params.Npulsars,gnx,gny,gnz); // Use when you want to run the simulation with the galactic potential + compute DM + geometry of emission
	cudaDeviceSynchronize();

	cudaMemcpy(x_bis,gx,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(y_bis,gy,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(z_bis,gz,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(gl_bis,ggl,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(gb_bis,ggb,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(vx_bis,gvx,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(vy_bis,gvy,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(vz_bis,gvz,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(error_bis,error_rel,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(DM_bis,gDM,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(w_r_fast_bis,gw_r_fast,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(w_r_pmps_bis,gw_r_pmps,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(xi_bis,gxi,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(rho_bis,grho,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(nx_bis,gnx,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(ny_bis,gny,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(nz_bis,gnz,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);

	cudaFree(gvx);cudaFree(gvy);cudaFree(gvz);cudaFree(gx);cudaFree(gy);cudaFree(gz);cudaFree(gage_pulsar);cudaFree(gsize);cudaFree(error_rel);cudaFree(gtmilky);cudaFree(ggb);cudaFree(ggl);cudaFree(gDM);cudaFree(gw_r_fast);cudaFree(gw_r_pmps);
	cudaFree(gperiod);cudaFree(galpha);cudaFree(g_nomegax);cudaFree(g_nomegay);cudaFree(g_nomegaz);cudaFree(gxi);cudaFree(grho);
	cudaFree(gnx);cudaFree(gny);cudaFree(gnz);

	for (int i=0;i<params.Npulsars;i++){
		params.x[i]=x_bis[i];
		params.y[i]=y_bis[i];
		params.z[i]=z_bis[i];
		params.vx[i]=vx_bis[i];
		params.vy[i]=vy_bis[i];
		params.vz[i]=vz_bis[i];
		params.err_rel_g[i]=error_bis[i];
		params.gl[i]=gl_bis[i];
		params.gb[i]=gb_bis[i];

                params.x_s[i]= params.x[i]; // shift center of the Galaxy to the Sun
                params.y_s[i]= params.y[i]-8.5;
                params.z_s[i]= params.z[i]-0.015;
                params.dist[i]= sqrt(sq(params.x_s[i])+sq(params.y_s[i])+sq(params.z_s[i]));
		params.DM[i]=DM_bis[i];
		params.w_r_fast[i]=w_r_fast_bis[i];
                params.w_r_pmps[i]=w_r_pmps_bis[i];
		params.xi[i]=xi_bis[i];
		params.rho[i]=rho_bis[i];
		params.nx[i]=nx_bis[i];
		params.ny[i]=ny_bis[i];
		params.nz[i]=nz_bis[i];
	}
	free(x_bis);free(y_bis);free(z_bis);free(vx_bis);free(vy_bis);free(vz_bis);free(error_bis);free(gl_bis);free(gb_bis);free(DM_bis);free(w_r_fast_bis);free(w_r_pmps_bis);free(xi_bis);free(rho_bis);
	free(nx_bis);free(ny_bis);free(nz_bis);
	cudaStatus = cudaGetLastError();
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "Erreur CUDA : %s\n", cudaGetErrorString(cudaStatus));
        // Gérer l'erreur
        }

        cudaDeviceSynchronize();
	
	//kick(&params);//birth kick velocity // Use when you don't want the galactic potential

	printf("##Total number of simulated pulsars: %ld \n",params.Npulsars);
	/* ISM modelling */
//        get_dispersion_measure(&params); 

	/* Computation of the number of orbits made by each NS  */
	nb_orbit(&params);

	/* Detection */
	//geometry(&params); // calculates the angles xi and w_r
        //pulse_profile_complete_2(&params); //Computes the pulse profile taking into account the DM + scattering + instrument
	sky_temp_Fmin_fermi(&params); //Get the sky temperatur at the position (l,b) of the pulsar. Maps of Haslam et al. (1982), reworked by Remazailles et al. (2015) + get the sensitivity of Fermi/LAT at a given position with python code from fermi 
	X_telescope_sky_coverage(&params); //allows to get info if the pulsar's position was observed by either XMM-Newton or Chandra
	radio_flux(&params); //calculates the radio flux of each pulsar
	get_fomega(&params); //Get all the different values of f_omega for the different angles of chi and zeta
	gamma_flux(&params); //idem for gamma flux
	X_flux(&params); //Idem Thermal X-ray flux
	check_x_pulse(&params); //Compute the X-ray pulsed fractions
	spinvel_angle(&params); //Computes the angle between the velocity vector and the rotation axis
	gamma_ray_peak_sep(&params); //Computes the gamma-ray peak separation 
	detection(&params); // check if the pulsar is beaming to us and if its flux is high enough to be detected
	detection_X(&params); //Idem for thermal X-ray flux
	//save_all(&params); //Save the info of every simulated pulsar


//Should I free all the parameters??         
	free(params.Pinit); // deallocate the buffer 
	free(params.period); // deallocate the buffer 
        free(params.Binit); // deallocate the buffer 
        free(params.alpha); // deallocate the buffer 
        free(params.x); // deallocate the buffer 
        free(params.y); // deallocate the buffer 
        free(params.z); // deallocate the buffer
	free(params.B);
        free(params.age_pulsar);
        free(params.gb);
        free(params.gl);
        free(params.x0);
        free(params.y0);
        free(params.z0);
        free(params.x_s);
        free(params.y_s);
        free(params.z_s);
        free(params.vx0);
        free(params.vy0);
        free(params.vz0);
        free(params.vx);
        free(params.vy);
        free(params.vz);
        free(params.err_rel_g);
        free(params.detec);
        free(params.detec_rad);
        free(params.detec_gam);
        free(params.detec_rg);
        free(params.Smin);
        free(params.Edot);
        free(params.Pdot);
        free(params.xi);
        free(params.rho);
        free(params.w_r);
        free(params.Fr);
        free(params.Fg);
        free(params.cos_a0);
        free(params.dist);
        free(params.n_omega_x);
        free(params.n_omega_y);
        free(params.n_omega_z);
	free(params.PA);
	free(params.DM);
	free(params.Nb_orb);
	free(params.delta);
	for(int i=0;i<91;i++) {free(params.fomega[i]);}
        free(params.fomega);
	free(params.Smin_fermi);
	free(params.w_int);
	free(params.w_r_fast);
        free(params.w_r_pmps);
	free(params.Smin_fast);
        free(params.Smin_pmps);
	free(params.temp);
	free(params.Fx);
	free(params.sky_XMM);
	free(params.sky_chandra);
	free(params.ex);
	free(params.ey);
	free(params.ez);
	free(params.n_mu_x);
	free(params.n_mu_y);
	free(params.n_mu_z);
	free(params.nx);
	free(params.ny);
	free(params.nz);
	free(params.cos_i);
	free(params.Temp);
	free(params.r_h);
	free(params.PF);

	gsl_rng_free(params.r);

	time(&end_time);
	int elapsed_time=difftime(end_time,start_time);
	double hour=((double)elapsed_time)/3600.0;double minu=(hour-((int)hour))*60;double sec=(minu-((int)minu))*60;
	int hour_s=(int)hour;int min_s=(int)minu;int sec_s=(int)sec;
	printf("Time taken for the simulation in h:m:s : %d:%d:%d\n",hour_s,min_s,sec_s);

return(0);


}

