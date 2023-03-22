#include"macro.h"
#include <stdlib.h>
#include<string.h>
#include <stdio.h>
#include"initialize.h"
#include<gsl/gsl_rng.h>
#include <math.h>
#include "birth_pulsars.h"
#include "evolution.h"
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_roots.h>



double func_angle_mhd(double x, void *params){ // function for which we want to find the roots Eq. 20 from Phillipov 2014
  
      struct func_params *part= (struct func_params*)params;
      double result;
      double k;
      long np=part->np;
      double alpha_d=part->alpha_d;
      double tau_d=part->tau_d;
      double age = part->age_pulsar[np];
      double cosa0 = part->cos_a0[np] ;
      double t_used;
      double sin2a0=1.-sq(cosa0);

	     if (part->Bfield_var==1) t_used  =   (alpha_d*tau_d)*(pow(1.+age/tau_d,1.-2/alpha_d)-1)/(alpha_d-2.); //Eq.19 from Dirson et al. 2022

	            else              t_used  =   part->age_pulsar[np];  // Eq.20 From Phillipov 2014

      k  = t_used/part->tau_MHD_al + 0.5/sin2a0 + log(sqrt(sin2a0));


      return result= 0.5 /sq(x) + log(x) - k;

}
/*

my_f(double x, void *params){



}

double my_df(double x, void *params)
{
  struct func_params *part= (struct func_params*)params;
  return -1./(x*x*x)+1./x; 

}

void my_fdf(double x, void *params, double *f,double *df){

  struct func_params *part= (struct func_params*)params;
  double cosa0 = part->cos_a0 ;
  double sin2a0 = 1.0-sq(cosa0) ;

  double k=part->age_pulsar[part->np]/part->tau_MHD_al + 0.5/sin2a0 + log(sqrt(sin2a0));
  double t  = 0.5 /sq(x) + log(x) - k;
  *f        = t;
  *df       = -1./(x*x*x)+1./x;


}

*/

/*****************************************************************  Temporal evolution of the pulsars *********************************************************************/

int evolution(void *params){  

        struct func_params *part= (struct func_params*)params;
  
        FILE *fp;
        if ( ( fp = fopen("pulsar_properties.dat","w+")) == NULL){
      	     printf("Couldn't open file pulsar_properties.dat \n");
	     exit(-1);
         }
        /*char *fname;
        fname=malloc(50);
        sprintf(fname,"population_init%.0ld.dat",part->Npulsars);
        */
        long np;
        long count=0;
        double k;
        double four_pi2 = 4*M_PI*M_PI;   
        double two_pi = 2*M_PI;   
        double height_pi2 =  8*M_PI*M_PI;   
        double k1 = 5.277e-33 ; // en 1/s (8*M_PI*pow(part->R,6))/(3*SI_mu0*cube(SI_C)*SI_I);
        double tau_0;
        double sina; 
        double alpha;
        double omega;
        double omega_0;
        double omega_dot;
  
  
        int status;
        int iter = 0;
        const gsl_root_fsolver_type *R;
        gsl_root_fsolver *s;
       //const gsl_root_fdfsolver_type *R;
      //gsl_root_fdfsolver *s;
  
  
	  for(np=0;np<part->Npulsars;np++){

		part->np = np; 
	        tau_0             =    1.0/k1 * sq(1.0/part->Binit[np])*sq(0.5*part->Pinit[np]/M_PI); //s-1
		part->cos_a0[np]  =    gsl_rng_uniform(part->r); // uniform between 0 and 1
		double cosa0      =    part->cos_a0[np] ;
		double cos2a0     =    cosa0*cosa0 ;
		omega_0           =    two_pi/(part->Pinit[np]);
		part->tau_MHD_al  =    tau_0*(1.0-cos2a0)/(cos2a0*cos2a0);
		part->tau_vac_al  =    1.5*tau_0 / cos2a0 ;
                part->tau_d       =    part->tau0_B0*pow(part->Binit[np],-part->alpha_d);

			if(part->ff_evol){  	//evolution of the angle alpha in the MHD case
  				double x_lo = 1.0e-17, x_hi = sqrt(1-cos2a0);
				//printf("f(x_low) %e f(x_high) %e \n",func_angle_mhd(x_lo,params),func_angle_mhd(x_hi,params));


				   if(fabs(func_angle_mhd(x_lo,params))<1e-15){

				   sina  = x_lo;	

				   } else if(fabs(func_angle_mhd(x_hi,params))<1e-15){

				   sina  = x_hi;	

				   } else{ 

				  
    
				        gsl_function  F;
				        F.function         =    &func_angle_mhd;
				        F.params           =    params;
				        R                  =    gsl_root_fsolver_brent;
				        s 		   =    gsl_root_fsolver_alloc(R);
					/* An other root-finding algorithm */ 
				        //gsl_function_fdf  FDF;
				        //FDF.f 	     =    &func_angle_mhd;
				        //FDF.df             =    &my_df;
				        //FDF.fdf            =    &my_fdf;
				        //FDF.params         =    params;
				        //R                    =    gsl_root_fdfsolver_newton;
				        //s 		      =    gsl_root_fdfsolver_alloc(R);

				      
				         gsl_root_fsolver_set(s,&F,x_lo,x_hi);      
				      
					      for(iter = 0; iter < 100; iter++) { // What happens if iter=100 is reached and no solutions found???
					
						    status = gsl_root_fsolver_iterate(s);
							
						    double left_int  = gsl_root_fsolver_x_lower(s);
						    double right_int = gsl_root_fsolver_x_upper(s);
							
					   	    //printf("iteration %03d: [%.010lf, %.010lf]\n", iter, left_int, right_int);
					
						    status = gsl_root_test_interval(left_int, right_int, 1.0e-10, 1.0e-15);
							if(status != GSL_CONTINUE) {
							  
							  //	    printf("status: %s\n", gsl_strerror(status));
				//			  	    printf("\nRoot interval = [%.010lf, %.010lf]\n", left_int, right_int);
							  sina=(left_int+right_int)/2.;
							  break;
							}
					      }

				          gsl_root_fsolver_free (s);

				    }

			    alpha            =   asin(sina);
			    part->alpha[np]  =   alpha;
			    omega            =   omega_0*(cos2a0*sina)/((sqrt(1-cos2a0))*sq(cos(alpha)));
			    part->period[np] = two_pi/omega;

			    if(part->Bfield_var==1){ 
			       part->B[np] = part->Binit[np]*pow(1+part->age_pulsar[np]/part->tau_d,-1./part->alpha_d);
			    } else {
			       part->B[np] = part->Binit[np];
			    }
	
			    omega_dot      = 1.5*k1*sq(part->B[np])*cube(omega)*(1+sq(sina)); // omega_dot = -K omega**3, K=(3/2)k1*B**2*(1+sina**2)
			    part->Pdot[np] = two_pi*omega_dot/sq(omega);
         	            part->Edot[np] =  SI_I*omega*omega_dot; 
			   
		       }

			
                                              
		    else if(part->vacuum){  // rewrite with omega and omega_dot
		      k   =   k1*(sq(part->Binit[np])*(1-sq(cosa0))); //vacuum
		      alpha = acos(cosa0);
	              part->alpha[np]=alpha;
		      part->period[np]  =   two_pi*sqrt(2.*(k*(part->age_pulsar[np])+sq(part->Pinit[np])/(height_pi2)));
	              part->Pdot[np]    =   four_pi2*k*(1./part->period[np]);
		    }
		    
		    else { 
		      k   =   k1*(sq(part->Binit[np])*((1-sq(cosa0))*exp(-2*part->age_pulsar[np]/part->tau_vac_al))); //vacuum
		      alpha = sqrt(1-sq(cosa0))*exp(-part->age_pulsar[np]/part->tau_vac_al);
	              part->alpha[np]=alpha;
		      part->period[np]=part->Pinit[np]*cos(alpha)/cosa0;
	              part->Pdot[np]    =   four_pi2*k*(1./part->period[np]);
		    }		
		    
	   


		      
	}	  

   fclose(fp);


return(0);

}
