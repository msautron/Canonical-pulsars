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



double func_angle_mhd(double x, void *params){
  
  struct func_params *part= (struct func_params*)params;
  double result;
  double k;
  long np=part->np;
  double alpha_d=part->alpha_d;
  double tau_d=part->tau_d;
  double cosa0 = part->cos_a0[part->np] ;
  double t_used;
  double sin2a0=1.-sq(cosa0);
// faire une condition où on choisit un t_utile
// utiliser tau_mhd

  if (part->Bfield_var==1) t_used=(alpha_d*tau_d)*(pow(1.+part->age_pulsar[np]/tau_d,1.-2/alpha_d)-1)/(alpha_d-2.);
	else  t_used=part->age_pulsar[np];
  k=t_used/part->tau_MHD_al + 0.5/sin2a0 + log(sqrt(sin2a0));
// k=part->age_pulsar[part->np]/part->tau_MHD_al + 0.5/sin2a0 + log(sqrt(sin2a0));
  //k=0.5001825;
  //printf("%e %e \n", part->age_pulsar[part->np], part->tau_MHD_al) ;
//  printf("k %e const_bvar %e term %e \n",k,const_bvar,0.5/sin2a0 + log(sqrt(sin2a0)));
//  printf("k= %.8f age/tau %e 1/(2sin2a0) %e log(sina0) %e \n",k,part->age_pulsar[part->np]/part->tau_MHD_al,0.5/sin2a0,log(sqrt(sin2a0)));
//  printf("k= %.8f age/tau %e age %e tau_MHD %e \n",k,part->age_pulsar[part->np]/part->tau_MHD_al,part->age_pulsar[part->np],part->tau_MHD_al);
 // printf("k= %.8f age/tau %e age %e tau_MHD %e \n",k,part->age_pulsar[part->np]/part->tau_MHD_al,part->age_pulsar[part->np],part->tau_MHD_al);
//  printf("cosa0 %.8f x %e \n",cosa0,x);
  result= 0.5 /sq(x) + log(x) - k;
  //printf("age_pulsar %e tau_mhd_al %e  cosa0 %e \n",part->age_pulsar[part->np],part->tau_MHD_al,part->cos_a0);
  return result;

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

int evolution(void *params){  
//	 gsl_rng * r;  /* global generator */
  struct func_params *part= (struct func_params*)params;
  //const gsl_rng_type * T;
  //       gsl_rng_env_setup();
  //     T = gsl_rng_default;
  //   r = gsl_rng_alloc (T);
  
  //FILE *fp;
  char *fname;
  fname=malloc(50);
  sprintf(fname,"population_init%.0ld.dat",part->Npulsars);
  //fp=fopen(fname,"r+");
  
  long np;
  long count=0;
  double k;
  double four_pi2=4*M_PI*M_PI;   
  double two_pi=2*M_PI;   
  double height_pi2=8*M_PI*M_PI;   
  double k1 = 5.277e-33 ; // en 1/s (8*M_PI*pow(part->R,6))/(3*SI_mu0*cube(SI_C)*SI_I);
  double tau_0;
  double sina; 
  double alpha;
  double omega;
  double omega_0;
  double omega_dot;
  //double tau0_B0 = 2*7e4*pow(3e8,part->alpha_d)*365*24*3600; //from Vigano
  //double p_age_sec; //age of the pulsar in s
  
  
  int status;
  int iter = 0;
  const gsl_root_fsolver_type *R;
  gsl_root_fsolver *s;
  //const gsl_root_fdfsolver_type *R;
  //gsl_root_fdfsolver *s;
  //double alpha;
  
 


 
	   //printf("Binit %e Binit %e \n",part->Binit[0],part->Pinit[0]);
  
  /* /!\ Problem when cosa0=0! This occurs for np=100 */
  
	  for(np=0;np<part->Npulsars;np++){
	//	if(fp!=NULL){
	//		fscanf(fp,"%lf %lf %lf %ld",&part->Binit[np],&part->Pinit[np],&part->age_pulsar[np],&np);
	//	 }
		part->np=np; 
	        tau_0  =   1.0/k1 * sq(1.0/part->Binit[np])*sq(0.5*part->Pinit[np]/M_PI); //s-1
		part->cos_a0[np]  =  gsl_rng_uniform(part->r); // uniform between 0 and 1
		double cosa0 = part->cos_a0[np] ;
		double cos2a0 = cosa0*cosa0 ;
		omega_0 = two_pi/(part->Pinit[np]);
		part->tau_MHD_al  =   tau_0*(1.0-cos2a0)/(cos2a0*cos2a0);
		part->tau_vac_al = 1.5*tau_0 / cos2a0 ;
                part->tau_d      = part->tau0_B0*pow(part->Binit[np],-part->alpha_d);
              //  part->tau_d      = part->tau0_B0*pow(2.5e8,-part->alpha_d)/365/24/3600;
		//part->tau_d = 1e9*365*24*3600;
	   //printf("tau_d %e \n",part->tau_d);
//		printf("%e %e \n",acos(cosa0),cosa0);	

//  printf("age/tau %e age %e tau_MHD %e \n",part->age_pulsar[part->np]/part->tau_MHD_al,part->age_pulsar[part->np],part->tau_MHD_al);
              //  printf("%e %e \n",part->tau_MHD_al/(365*24*3600),cosa0);
			if(part->ff_evol){  	//evolution of the angle alpha in the MHD case
  				double x_lo = 1.0e-17, x_hi = sqrt(1-cos2a0);
				//printf("f(x_low) %e f(x_high) %e \n",func_angle_mhd(x_lo,params),func_angle_mhd(x_hi,params));
            		        //printf("alpha0 %.9f \n",acos(cosa0));


				   if(fabs(func_angle_mhd(x_lo,params))<1e-15){

				   sina  = x_lo;	
                          //         printf("cas 1");

				   } else if(fabs(func_angle_mhd(x_hi,params))<1e-15){

				   sina  = x_hi;	
                        //           printf("cas 2");

				   } else{ 

				  
				    // printf("t/tau %e 1/2sin^2(x) %e log(sin(x)) %e \n",part->age_pulsar/part->tau_MHD_al,1./(2*(1-sq(part->cos_a0))),log(sqrt(1-sq(part->cos_a0))));    
    
				      gsl_function  F;
				      F.function         =    &func_angle_mhd;
				      F.params           =    params;
				      //gsl_function_fdf  FDF;
				      //FDF.f 	         =    &func_angle_mhd;
				      //FDF.df             =    &my_df;
				      //FDF.fdf            =    &my_fdf;
				      //FDF.params           =    params;
				      R                  =    gsl_root_fsolver_brent;
				      s 		 =    gsl_root_fsolver_alloc(R);
				      //R                  =    gsl_root_fdfsolver_newton;
				      //s 		 =    gsl_root_fdfsolver_alloc(R);

				      //	printf("part->cos_a0 %e \n",part->cos_a0);
				      //printf("%e %e \n", tau_0, part->tau_MHD_al) ;
				      
				      gsl_root_fsolver_set(s,&F,x_lo,x_hi);      
				      
					      for(iter = 0; iter < 100; iter++) { // QUE SE PASSE-T-IL SI ON ARRIVE a 100 ITERATIONS ET PAS DE SOLUTION ?
					
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

				    alpha = asin(sina);
				    part->alpha[np]=alpha;
				    omega = omega_0*(cos2a0*sina)/((sqrt(1-cos2a0))*sq(cos(alpha)));
				    part->period[np] = two_pi/omega;
				    //part->period[np]=(part->Pinit[np]*(sqrt(1-cos2a0))*sq(cos(alpha)))/(cos2a0*sina); 

				    if(part->Bfield_var==1){ 
			   	       part->B[np] = part->Binit[np]*pow(1+part->age_pulsar[np]/part->tau_d,-1./part->alpha_d);
			            } else {
			               part->B[np] = part->Binit[np];
				    }
		
				    //k     =   1.5*k1*(sq(part->Binit[np])*(1+sq(sina))); //force-free
				    omega_dot = 1.5*k1*sq(part->B[np])*cube(omega)*(1+sq(sina)); // omega_dot = -K omega**3, K=(3/2)k1*B**2*(1+sina**2)
				    part->Pdot[np] = two_pi*omega_dot/sq(omega);
				    if (part->period[np]<0.){
					printf("P is negative!!! Pinit %e P %e np %ld omega %e omega0 %e cos2a0 %e sina %e cosa %e sin2a0 %e \n",part->Pinit[np],part->period[np],np,omega,omega_0,cos2a0,sina,cos(alpha),sqrt(1-cos2a0));
				//	break;
					}
	                         //   part->Pdot[np]    =   four_pi2*k*(1./part->period[np]);
			//	printf("Pinf %e Pinit %e P %e Pdot %e tau_d %e \n",sqrt(square_Pinf),part->Pinit[np],part->period[np],part->Pdot[np],part->tau_d/(365*24*3600));	     
				//printf("Pinf %e Pinit %e P %e \n",sqrt(square_Pinf),part->Pinit[np],part->period[np]);	     
				//printf("tau_d %e B0 %e B %e \n",part->tau_d/(365*24*3600),part->Binit[np],part->B[np]);	     
			//	printf("tau0B0 %e Binit %e B0^alpha %e \n",tau0_B0,part->Binit[np],pow(part->Binit[np],-alpha_d));	     
            //printf("np %ld cosa0 %.9f cos(alpha) %.9f age %.9f \n",np,part->cos_a0[np],cos(part->alpha[np]),part->age_pulsar[np]/(365*24*3600));
            //printf("%.9f %.9f \n",part->cos_a0[np],cos(part->alpha[np]));
			   
		  }
    
		    else if(part->vacuum){
		      k   =   k1*(sq(part->Binit[np])*(1-sq(cosa0))); //vacuum
		      alpha = acos(cosa0);
	              part->alpha[np]=alpha;
		      part->period[np]  =   two_pi*sqrt(2.*(k*(part->age_pulsar[np])+sq(part->Pinit[np])/(height_pi2)));
	              part->Pdot[np]    =   four_pi2*k*(1./part->period[np]);
		    }
		    
		    else if(part->vacuum_evol){ 
		      k   =   k1*(sq(part->Binit[np])*((1-sq(cosa0))*exp(-2*part->age_pulsar[np]/part->tau_vac_al))); //vacuum
		      alpha = sqrt(1-sq(cosa0))*exp(-part->age_pulsar[np]/part->tau_vac_al);
	              part->alpha[np]=alpha;
		      part->period[np]=part->Pinit[np]*cos(alpha)/cosa0;
	              part->Pdot[np]    =   four_pi2*k*(1./part->period[np]);
		    }		
		    
		    else {
		      k   =   (3/2.)*k1*(sq(part->Binit[np])*(1+1-sq(cosa0))); //force-free
		      alpha = acos(cosa0);
	              part->alpha[np]=alpha;
		      //printf("alpha %e \n",part->alpha[np]);
		    }  
	   
	    //part->alpha[np] = acos(cosa0);


//	    part->Edot[np]    =   (four_pi2*SI_I*part->Pdot[np])/cube(part->period[np]); // in J/s
	    part->Edot[np]    =  SI_I*omega*omega_dot; 
	   //printf("np %e age %.10e P %e Pinit %e Binit %e Pdot %e Edot %e \n",(double)np,part->age_pulsar[np],part->period[np],part->Pinit[np],part->Binit[np],part->Pdot[np],part->Edot[np]);
		      
	 //printf(" %e %e \n",cosa0,part->tau_MHD_al/(365*3600*24));
	 //	printf("a %e a0 %e tau_vac_al %e \n",part->alpha[np],acos(part->cos_a0),part->tau_vac_al);
            //printf("%ld %.9f %.9f %.9f \n",np,part->cos_a0[np],cos(part->alpha[np]),part->age_pulsar[np]/(365*24*3600));
            //printf("%.9f %.9f \n",cos(part->alpha[np]),part->age_pulsar[np]/(365*24*3600));
		      	//printf(" %e  %e \n",part->period[np],part->Pdot[np]);
					//printf("%e %e %e %e \n",acos(part->cos_a0[np])*180/M_PI,part->alpha[np]*180/M_PI,part->cos_a0[np],cos(part->alpha[np]));
//printf("np %e age %.10e P %e Pinit %e Binit %e Pdot %e Edot %e B %e \n",(double)np,part->age_pulsar[np],part->period[np],part->Pinit[np],part->Binit[np],part->Pdot[np],part->Edot[np],part->B[np]);
//printf("P %e Pdot %e \n",part->period[np],part->Pdot[np]);

	    //printf("np %e age %.10e P %e Pinit %e Binit %e Pdot %e Edot %e \n",(double)np,part->age_pulsar[np],part->period[np],part->Pinit[np],part->Binit[np],part->Pdot[np],part->Edot[np]);
		    if(part->Edot[np]>1e31){  //1e28 J/s is equal to 1e35 erg/s

//			      printf("alpha0 %.9f alpha %.9f numerateur %e denominateur %e \n",acos(cosa0), alpha,(part->Pinit[np]*(sqrt(1-cos2a0))*sq(cosa)),(cos2a0*sina));
//			      printf("cosa0 %.9f cos2a0 %.9f sina0 %.9f sina %.9f \n",cosa0,cos2a0,sqrt(1-cos2a0),sin(part->alpha[np]));
		      //        would be better to print it in a file
			      //             printf(" %e \n",k);
			      count++;
			      //	printf("p_age_sec %e tau_0 %e tau_vac_al %e k %e \n",p_age_sec,tau_0,tau_vac_al,k);
		    }	
	}	  
//	printf("Number of pulsars with Edot > 1e35 erg.s-1: %ld \n",count);

//fclose(fp);
//gsl_rng_free (part->r);


return(0);

}
