#define GNU_SOURCE
//#define _GNU_SOURCE
#include <stdlib.h>
#include"initialize.h"
#include"macro.h"
#include"birth_pulsars.h"
#include"evolution.h"
#include"galac_pot.h"
#include"ism_scattering.h"
#include"detection.h"
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#define TAILLE_MAX 1000 // Tableau de taille 1000

#define NTHREAD (500)
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

__global__ void evol_galac_PEFRL(double *gx,double *gy,double *gz, double *gvx,double *gvy,double *gvz,double *gage_pulsar,int *gsize,double *error_rel,double *gtmilky){ //PEFRL integration scheme

    const double yr_sec=365*24*3600;int np;
    double step;
    const double kpc2km=3.0856775807e16;
    double L0=1.0;double v0=100.0;
    double T0=(L0*kpc2km)/v0;
    double step_wd;
    double xsi=0.1786178958448091;double lambda=-0.2123418310626054;
    double chi=-0.6626458266981849e-01;
    int idx=threadIdx.x+blockIdx.x*NTHREAD;
    int size=*gsize;
    int start=idx*size;int end=start+size;
    int count=0;
    //printf("Bonjour I'm thread %d in block %d, with absolute %d\n",threadIdx.x,blockIdx.x,idx);
    //for(np=start;np<end;np++){
    while(count<size){

       printf("Count=%d\n",count);
       if (count%2==0) np=(long)((idx*size+count)/2);
       else if (count%2==1) np=(long)(10000000-((idx*size+count+1)/2));

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

      count+=1;
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

	printf("## birth_rate (yr) %ld \n",params.birth_rate);
	printf("## sigma_P (s or log(s)) %e \n",params.sigma_p);
	printf("## sigma_B (logB) %e \n",params.sigma_b);
	printf("## Pmean (s) %e \n",params.p_mean);
	printf("## Bmean (T) %e \n",params.b_mean);
	printf("## k_tau0_B0  %ld \n",params.k_tau0_B0);
	printf("## alpha_d  %e \n",params.alpha_d);
	printf("## vold (km.s-1) %ld \n",params.v_old);
	printf("## vyoung (km.s-1) %ld \n",params.v_young);
	printf("## B var %ld \n",params.Bfield_var);
	printf("\n");

         /* Calculates P, dotP and B */
	birth(&params); // generates the pulsar with an initial Pinit, Binit and age
        evolution(&params); //calculates P,Pdot and Edot
  
        // Distribution of the pulsars in the Galaxy
	distrib_init(&params);    // initial distribution of the pulsars    
	//distrib_init_2(&params); // initial distribution of the pulsars, when they are born in spiral arms
	distrib_vinit(&params); // Use when you want to run the simulation with the galactic potential
	
	double *gvx;double *gvy;double *gvz;double *gx;double *gy;double *gz;double *gage_pulsar;int *gsize;double *error_rel;double *gtmilky;
	double *x_bis;double *y_bis;double *z_bis;double *vx_bis;double *vy_bis;double *vz_bis;double *error_bis;double tmilky;
	const double yr_sec=365*24*3600;tmilky=TMILKY*yr_sec;
	int size=params.Npulsars/(NBLOCK*NTHREAD);
        cudaMalloc((void**)&gvx,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gvy,sizeof(double)*params.Npulsars);
 	cudaMalloc((void**)&gvz,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gx,sizeof(double)*params.Npulsars);
 	cudaMalloc((void**)&gy,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gz,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gage_pulsar,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gsize,sizeof(int));
	cudaMalloc((void**)&error_rel,sizeof(double)*params.Npulsars);
	cudaMalloc((void**)&gtmilky,sizeof(double));

	x_bis=(double *)calloc(params.Npulsars,sizeof(double));
	y_bis=(double *)calloc(params.Npulsars,sizeof(double));
	z_bis=(double *)calloc(params.Npulsars,sizeof(double));
	vx_bis=(double *)calloc(params.Npulsars,sizeof(double));
	vy_bis=(double *)calloc(params.Npulsars,sizeof(double));
	vz_bis=(double *)calloc(params.Npulsars,sizeof(double));
	error_bis=(double *)calloc(params.Npulsars,sizeof(double));

        cudaMemcpy(gvx,params.vx,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gvy,params.vy,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gvz,params.vz,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gx,params.x,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gy,params.y,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gz,params.z,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gage_pulsar,params.age_pulsar,params.Npulsars*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(gsize,&size,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(gtmilky,&tmilky,sizeof(double),cudaMemcpyHostToDevice);

        evol_galac_PEFRL<<<dimgrid,dimblock>>>(gx,gy,gz,gvx,gvy,gvz,gage_pulsar,gsize,error_rel,gtmilky); // Use when you want to run the simulation with the galactic potential
	cudaDeviceSynchronize();

	cudaMemcpy(x_bis,gx,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(y_bis,gy,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(z_bis,gz,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(vx_bis,gvx,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(vy_bis,gvy,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(vz_bis,gvz,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(error_bis,error_rel,params.Npulsars*sizeof(double),cudaMemcpyDeviceToHost);

	cudaFree(gvx);cudaFree(gvy);cudaFree(gvz);cudaFree(gx);cudaFree(gy);cudaFree(gz);cudaFree(gage_pulsar);cudaFree(gsize);cudaFree(error_rel);cudaFree(gtmilky);

	for (int i=0;i<params.Npulsars;i++){
		params.x[i]=x_bis[i];
		params.y[i]=y_bis[i];
		params.z[i]=z_bis[i];
		params.vx[i]=vx_bis[i];
		params.vy[i]=vy_bis[i];
		params.vz[i]=vz_bis[i];
		params.err_rel_g[i]=error_bis[i];

                params.x_s[i]= params.x[i]; // shift center of the Galaxy to the Sun
                params.y_s[i]= params.y[i]-8.5;
                params.z_s[i]= params.z[i]-0.015;
                params.dist[i]= sqrt(sq(params.x_s[i])+sq(params.y_s[i])+sq(params.z_s[i]));
                double r= sqrt(sq(params.x_s[i])+sq(params.y_s[i]));
	        double gbr;double glr;double RAD=180/M_PI;
                if(params.dist[i]<1e-15){
                   gbr=0;
                }else{
                   gbr=asin(params.z_s[i]/params.dist[i]);
                }
                params.gb[i]=gbr*RAD;

                if(r<1e-15){
                   glr=0;
                  }
                else{
                    if(params.x_s[i]>=0){
                         glr=acos(-params.y_s[i]/r);
                                 }
                    else{
                         glr=acos(params.y_s[i]/r)+M_PI;
                    }
                }

                params.gl[i]=glr*RAD;
	}
	free(x_bis);free(y_bis);free(z_bis);free(vx_bis);free(vy_bis);free(vz_bis);free(error_bis);
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

	/* Detection */
	geometry(&params); // calculates the angles xi and w_r
	radio_flux(&params); //calculates the radio flux of each pulsar
	gamma_flux(&params); //idem for gamma flux
	radio_flux_low_freq(&params);
	spinvel_angle(&params);
	detection(&params); // check if the pulsar is beaming to us and if its flux is high enough to be detected



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
        free(params.flux_low_freq);
        free(params.Fg);
        free(params.cos_a0);
        free(params.dist);
        free(params.n_omega_x);
        free(params.n_omega_y);
        free(params.n_omega_z);
	free(params.PA);

	gsl_rng_free(params.r);

	time(&end_time);
	int elapsed_time=difftime(end_time,start_time);
	double hour=((double)elapsed_time)/3600.0;double minu=(hour-((int)hour))*60;double sec=(minu-((int)minu))*60;
	int hour_s=(int)hour;int min_s=(int)minu;int sec_s=(int)sec;
	printf("Time taken for the simulation in h:m:s : %d:%d:%d\n",hour_s,min_s,sec_s);

return(0);


}

