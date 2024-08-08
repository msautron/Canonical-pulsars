#include"macro.h"
#include <stdlib.h>
#include<string.h>
#include <stdio.h>
#include"initialize.h"
#include<gsl/gsl_rng.h>
#include <math.h>
#include "birth_pulsars.h"
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>


double f_weibull_P0(double P0,double J,double K){
	return (J/(K*sq(P0)))*(1.0/pow(K*P0,J-1))*exp(-pow(1.0/(K*P0),J));
}

int birth(void *params){//generates Npulsars with an initial period, B and age 

        struct func_params *part= (struct func_params*)params;
	bool sample=false;
	double J=1.2;
	double K=14.5;
	double P0;
	double pdf_val,comp_val;
        long np=0;

           	while(np<part->Npulsars){  

                        //part->Pinit[np]=part->p_mean+gsl_ran_gaussian_ziggurat(part->r,part->sigma_p); //Normal distribution from most of the litterature
			part->Pinit[np]=pow(10,log10(part->p_mean)+gsl_ran_gaussian_ziggurat(part->r,part->sigma_p)); //Log normal distribution found in Igoshev et al. (2022)
			//if (part->Pinit[np] < 0){
		        //    continue;	
			//}

		        //Weibull distribution for the spin period Du et al. (2024) 
		        /*while (sample==false){
                                P0=(200e-3-5e-3)*gsl_rng_uniform(part->r)+5e-3;
                                pdf_val=f_weibull_P0(P0,J,K);
                                comp_val=f_weibull_P0(5e-3,J,K);
                                if (comp_val<=pdf_val) {sample=true;}
                        }
                        sample=false;
                        part->Pinit[np]=P0; */

			//Magnetic field initialization
                        part->Binit[np]=pow(10,log10(part->b_mean)+gsl_ran_gaussian_ziggurat(part->r,part->sigma_b));
			//Age initilization
		        part->age_pulsar[np]  =   part->birth_rate*np*365*24*3600; //s          
		        np++;
                }
       return(0);

        }
           
