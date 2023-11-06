#define GNU_SOURCE
#define _GNU_SOURCE
#define _GNU_SOURCE
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
#define TAILLE_MAX 1000 // Tableau de taille 1000

/*double period_evol(double t, void *params){
 *         double result;
 *                 struct func_params *part = (struct func_params *)params;
 *                         result=sqrt(-2/(3*cube(SI_C))*(cube(part->R)*part->)*t-125)
 *
 *                         }
 *                         */
int main(int argc, char **argv){

        struct func_params params;
    	params.T = gsl_rng_default;
    	params.r = gsl_rng_alloc (params.T);



	 time_t result = time(NULL);
         if(result != (time_t)(-1)) gsl_rng_set(params.r,(intmax_t)result);
             printf("# The current time is %s(%jd seconds since the Epoch)\n",
                  asctime(gmtime(&result)), (intmax_t)result);

	     //gsl_rng_set(params.r,123);  

        initialize(argc, argv, &params); //parameters initialization


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
	distrib_init(&params); // initial distribution of the pulsars
	//distrib_vinit(&params); // Use when you want to run the simulation with the galactic potential
        //evol_galac_PEFRL(&params); // Use when you want to run the simulation with the galactic potential
	kick(&params);//birth kick velocity // Use when you don't want the galactic potential

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


return(0);


}

