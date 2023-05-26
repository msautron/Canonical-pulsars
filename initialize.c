#include"initialize.h"
#include"macro.h"
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
void initialize(int argc, char *argv[], void *params){
        struct func_params *part= (struct func_params *)params;


	 if (argc>1){
                part->Npulsars   = (long) strtod(argv[1],(char **)NULL);
                part->sigma_b    = (double) strtod(argv[2],(char **)NULL);
                part->b_mean     = (double) strtod(argv[3],(char **)NULL);
                part->p_mean     = (double) strtod(argv[4],(char **)NULL);
                part->sigma_p    = (double) strtod(argv[5],(char **)NULL);
                part->alpha_d    = (double) strtod(argv[6],(char **)NULL);
                part->k_tau0_B0  = (double) strtod(argv[7],(char **)NULL);
                part->birth_rate = (double) strtod(argv[8],(char **)NULL);
                part->v_old      = (double) strtod(argv[9],(char **)NULL);
		//printf("sigma_b %e \n",part->sigma_b);
//                part->n_init= (long) strtod(argv[2],(char **)NULL);
         }
         else{
         part->Npulsars         =  	10000000;
         //part->Npulsars         =  	10000;
	 part->k_tau0_B0        =       5;
         part->birth_rate	= 	70;
         part->b_mean		= 	2.5e8; //Tesla usual value used 
	 //part->b_mean           =       275422870.33381635; //Tesla, value used in Igoshev et al. (2022)
         part->p_mean		= 	60e-3;// usual value used in seconds (normal distribution)
	 //part->p_mean           =       0.09120108393559097;// value used in Igoshev et al. (2022) in seconds (log normal distribution)
         part->sigma_p		= 	0.010;//usual value used in s (normal distribution)
	 //part->sigma_p          =       0.53; // value used in Igoshev et al. (2022) (log normal distribution)
         part->alpha_d          = 	1.5;
         //part->sigma_b		= 	0.5;//usual value used
	 part->sigma_b          =       0.5; //value used in Igoshev et al. (2022)
         part->v_old            =	75.; //  km/s 
/*
         part->birth_rate	= 	300;
         part->b_mean		= 	1.2e8; //Tesla
         part->p_mean		= 	50e-3;//s 
         part->sigma_p		= 	0.001;//s 
         part->sigma_b		= 	0.4 ;// 
         part->v_old            =	75.; //  km/s 
*/
  //       part->n_init=0.;
         }

     
       part->tau0_B0            = 	part->k_tau0_B0*7e4*pow(3e8,part->alpha_d)*365*24*3600; //from Vigano
       part->R			=	12000;//m
       part->zexp               =	0.5; // kpc
       part->Rexp		=	4.5; //kpc
       part->sigma_v		=	60.; //  km/s 
       part->v_young		=	265.; //  km/s 
       part->vacuum		=	0;
       part->vacuum_evol	=	0;
       part->ff_evol            =       1;
       part->NenuFAR            =       0;
       part->ska                =       0;
       part->Bfield_var         =       1;


       part->Pinit = malloc(sizeof(*part->Pinit) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Pinit == NULL) printf("Pinit: allocation failed"); // check if allocation succeeded 
       part->Binit = malloc(sizeof(*part->Binit) * part->Npulsars); // initialize pointer (allocate) 
                if (part->Binit == NULL) printf("Binit: allocation failed"); // check if allocation succeeded 
       part->alpha= malloc(sizeof(*part->alpha) * part->Npulsars); // initialize pointer (allocate) 
                if (part->alpha== NULL) printf("alpha: allocation failed"); // check if allocation succeeded 
       part->age_pulsar= malloc(sizeof(*part->age_pulsar) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->age_pulsar == NULL) printf("t_pulsar: allocation failed"); // check if allocation succeeded 
       part->period= malloc(sizeof(*part->period) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->period == NULL) printf("period: allocation failed"); // check if allocation succeeded 
       part->Edot= malloc(sizeof(*part->Edot) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Edot== NULL) printf("Edot: allocation failed"); // check if allocation succeeded 
       part->dist= malloc(sizeof(*part->dist) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->dist== NULL) printf("dist: allocation failed"); // check if allocation succeeded 
       part->x= malloc(sizeof(*part->x) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x== NULL) printf("x: allocation failed"); // check if allocation succeeded 
       part->y= malloc(sizeof(*part->y) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y== NULL) printf("y: allocation failed"); // check if allocation succeeded 
       part->z= malloc(sizeof(*part->z) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z== NULL) printf("z: allocation failed"); // check if allocation succeeded 
       part->Pdot= malloc(sizeof(*part->Pdot) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Pdot== NULL) printf("Pdot: allocation failed"); // check if allocation succeeded 
       part->Fr= malloc(sizeof(*part->Fr) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fr== NULL) printf("Fr: allocation failed"); // check if allocation succeeded 
       part->Fg= malloc(sizeof(*part->Fg) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fg== NULL) printf("Fg: allocation failed"); // check if allocation succeeded 

       part->x0= malloc(sizeof(*part->x0) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x0== NULL) printf("x0: allocation failed"); // check if allocation succeeded 
       part->y0= malloc(sizeof(*part->y0) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y0== NULL) printf("y0: allocation failed"); // check if allocation succeeded 
       part->z0= malloc(sizeof(*part->z0) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z0== NULL) printf("z0: allocation failed"); // check if allocation succeeded 
       part->xi= malloc(sizeof(*part->xi) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->xi== NULL) printf("xi: allocation failed"); // check if allocation succeeded 
       part->rho= malloc(sizeof(*part->rho) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->rho== NULL) printf("rho: allocation failed"); // check if allocation succeeded 
       part->w_r= malloc(sizeof(*part->w_r) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_r== NULL) printf("w_r: allocation failed"); // check if allocation succeeded 
       part->Smin= malloc(sizeof(*part->Smin) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin== NULL) printf("Smin: allocation failed"); // check if allocation succeeded 
       part->cos_a0= malloc(sizeof(*part->cos_a0) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->cos_a0== NULL) printf("cos_a0: allocation failed"); // check if allocation succeeded 
       part->B= malloc(sizeof(*part->B) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->B== NULL) printf("B: allocation failed"); // check if allocation succeeded 
       part->flux_low_freq= malloc(sizeof(*part->flux_low_freq) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->flux_low_freq== NULL) printf("flux_low_freq: allocation failed"); // check if allocation succeeded 
       part->gl= malloc(sizeof(*part->gl) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->gl== NULL) printf("gl: allocation failed"); // check if allocation succeeded 
       part->gb= malloc(sizeof(*part->gb) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->gb== NULL) printf("gb: allocation failed"); // check if allocation succeeded 

       part->vx0= malloc(sizeof(*part->vx0) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->vx0== NULL) printf("vx0: allocation failed");	// check if allocation succeede

       part->vy0= malloc(sizeof(*part->vy0) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vy0== NULL) printf("vy0: allocation failed"); // check if allocation succeede

       part->vz0= malloc(sizeof(*part->vz0) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vz0== NULL) printf("vz0: allocation failed"); // check if allocation succeede
       part->vz= malloc(sizeof(*part->vz) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vz== NULL) printf("vz: allocation failed"); // check if allocation succeede

       part->vy= malloc(sizeof(*part->vy) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vy== NULL) printf("vy: allocation failed"); // check if allocation succeede
       part->vx= malloc(sizeof(*part->vx) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vx== NULL) printf("vx: allocation failed"); // check if allocation succeede

       part->err_rel_g= malloc(sizeof(*part->err_rel_g) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->err_rel_g== NULL) printf("err_rel_g: allocation failed"); // check if allocation succeede
       part->x_s= malloc(sizeof(*part->x_s) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x_s== NULL) printf("x_s: allocation failed"); // check if allocation succeeded
       part->y_s= malloc(sizeof(*part->y_s) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y_s== NULL) printf("y_s: allocation failed"); // check if allocation succeeded  
       part->z_s= malloc(sizeof(*part->z_s) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z_s== NULL) printf("z_s: allocation failed"); // check if allocation succeeded 
       part->detec= malloc(sizeof(*part->detec) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec== NULL) printf("detec: allocation failed"); // check if allocation succeeded
       part->detec_rad= malloc(sizeof(*part->detec_rad) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_rad== NULL) printf("detec_rad: allocation failed"); // check if allocation succeeded
       part->detec_gam= malloc(sizeof(*part->detec_gam) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_gam== NULL) printf("detec_gam: allocation failed"); // check if allocation succeeded
       part->detec_rg= malloc(sizeof(*part->detec_rg) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_rg== NULL) printf("detec_rg: allocation failed"); // check if allocation succeeded
       part->n_omega_z= malloc(sizeof(*part->n_omega_z) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_z== NULL) printf("n_omega_z: allocation failed"); // check if allocation succeed
       part->n_omega_y= malloc(sizeof(*part->n_omega_y) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_y== NULL) printf("n_omega_y: allocation failed"); // check if allocation succeeded
       part->n_omega_x= malloc(sizeof(*part->n_omega_x) * part->Npulsars); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_x== NULL) printf("n_omega_x: allocation failed"); // check if allocation succeeded

}
