#include"initialize.h"
#include"macro.h"
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
void initialize(int argc, char *argv[], void *params){
        struct func_params *part= (struct func_params *)params;


	 if (argc>1){
                part->sigma_b    = (double) strtod(argv[1],(char **)NULL);
                part->b_mean     = (double) strtod(argv[2],(char **)NULL);
                part->p_mean     = (double) strtod(argv[3],(char **)NULL);
                part->sigma_p    = (double) strtod(argv[4],(char **)NULL);
                part->birth_rate = (double) strtod(argv[5],(char **)NULL);
		part->A_propto = (double) strtod(argv[6],(char **)NULL);
		part->D_propto = (double) strtod(argv[7],(char **)NULL);
		part->M_for_K = (double) strtod(argv[8],(char **)NULL);
		part->R_for_K = (double) strtod(argv[9],(char **)NULL);
		part->tau_d = (double) strtod(argv[10],(char **)NULL);
		part->alpha_d = (double) strtod(argv[11],(char **)NULL);
         }
         else{
         part->M_for_K                 =       1.4;
	 part->R_for_K                 =       R_NS;
	 part->A_propto                =       4.0e8;
	 part->D_propto                =       3.8e4;
         part->birth_rate	= 	34;
	 part->b_mean           =       275422870.33381635; //Tesla, value used in Igoshev et al. (2022)
	 part->p_mean           =       129e-3;//1.174898e-1;// value used in Igoshev et al. (2022) in seconds (log normal distribution)
	 part->sigma_p          =       0.45; // value used in Igoshev et al. (2022) (log normal distribution)
	 part->sigma_b          =       0.5; //value used in Igoshev et al. (2022)
	 part->alpha_d          =       1.5;
	 part->tau_d            =       1.8e6*365*24*3600; // in s 
         }

       part->Npulsars           =       1.5e6;//9e7/part->birth_rate;
       part->v_old              =       265.0;//km/s
       part->R			=	12000;//m
       part->zexp               =	0.18; // kpc
       part->Rexp		=	4.5; //kpc
       part->sigma_v		=	265.; //  km/s 
       part->v_young		=	265.; //  km/s 
       part->vacuum		=	0;
       part->vacuum_evol	=	0;
       part->ff_evol            =       1;
       part->NenuFAR            =       0;
       part->ska                =       0;
       part->Bfield_var         =       1;


       part->Pinit =(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Pinit == NULL) printf("Pinit: allocation failed"); // check if allocation succeeded 
       part->Binit = (double *)calloc(part->Npulsars,sizeof(double)); // initialize pointer (allocate) 
                if (part->Binit == NULL) printf("Binit: allocation failed"); // check if allocation succeeded 
       part->alpha= (double *)calloc(part->Npulsars,sizeof(double)); // initialize pointer (allocate) 
                if (part->alpha== NULL) printf("alpha: allocation failed"); // check if allocation succeeded 
       part->age_pulsar= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->age_pulsar == NULL) printf("t_pulsar: allocation failed"); // check if allocation succeeded 
       part->period= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->period == NULL) printf("period: allocation failed"); // check if allocation succeeded 
       part->Edot= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Edot== NULL) printf("Edot: allocation failed"); // check if allocation succeeded 
       part->dist= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->dist== NULL) printf("dist: allocation failed"); // check if allocation succeeded 
       part->x= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x== NULL) printf("x: allocation failed"); // check if allocation succeeded 
       part->y= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y== NULL) printf("y: allocation failed"); // check if allocation succeeded 
       part->z= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z== NULL) printf("z: allocation failed"); // check if allocation succeeded 
       part->Pdot= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Pdot== NULL) printf("Pdot: allocation failed"); // check if allocation succeeded 
       part->Fr= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fr== NULL) printf("Fr: allocation failed"); // check if allocation succeeded 
       part->Fg= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fg== NULL) printf("Fg: allocation failed"); // check if allocation succeeded 

       part->x0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x0== NULL) printf("x0: allocation failed"); // check if allocation succeeded 
       part->y0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y0== NULL) printf("y0: allocation failed"); // check if allocation succeeded 
       part->z0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z0== NULL) printf("z0: allocation failed"); // check if allocation succeeded 
       part->xi= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->xi== NULL) printf("xi: allocation failed"); // check if allocation succeeded 
       part->rho= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->rho== NULL) printf("rho: allocation failed"); // check if allocation succeeded 
       part->w_r= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_r== NULL) printf("w_r: allocation failed"); // check if allocation succeeded 
       part->Smin= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin== NULL) printf("Smin: allocation failed"); // check if allocation succeeded 
       part->cos_a0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->cos_a0== NULL) printf("cos_a0: allocation failed"); // check if allocation succeeded 
       part->B= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->B== NULL) printf("B: allocation failed"); // check if allocation succeeded 
       part->gl= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->gl== NULL) printf("gl: allocation failed"); // check if allocation succeeded 
       part->gb= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->gb== NULL) printf("gb: allocation failed"); // check if allocation succeeded 

       part->vx0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->vx0== NULL) printf("vx0: allocation failed");	// check if allocation succeede

       part->vy0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vy0== NULL) printf("vy0: allocation failed"); // check if allocation succeede

       part->vz0= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vz0== NULL) printf("vz0: allocation failed"); // check if allocation succeede
       part->vz= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vz== NULL) printf("vz: allocation failed"); // check if allocation succeede

       part->vy= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vy== NULL) printf("vy: allocation failed"); // check if allocation succeede
       part->vx= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->vx== NULL) printf("vx: allocation failed"); // check if allocation succeede

       part->err_rel_g= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->err_rel_g== NULL) printf("err_rel_g: allocation failed"); // check if allocation succeede
       part->x_s= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->x_s== NULL) printf("x_s: allocation failed"); // check if allocation succeeded
       part->y_s= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->y_s== NULL) printf("y_s: allocation failed"); // check if allocation succeeded  
       part->z_s= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->z_s== NULL) printf("z_s: allocation failed"); // check if allocation succeeded 
       part->detec= (long *)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec== NULL) printf("detec: allocation failed"); // check if allocation succeeded
       part->detec_rad=(long *)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_rad== NULL) printf("detec_rad: allocation failed"); // check if allocation succeeded
       part->detec_gam= (long*)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_gam== NULL) printf("detec_gam: allocation failed"); // check if allocation succeeded
       part->detec_rg= (long *)calloc(part->Npulsars,sizeof(long)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->detec_rg== NULL) printf("detec_rg: allocation failed"); // check if allocation succeeded
       part->n_omega_z= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_z== NULL) printf("n_omega_z: allocation failed"); // check if allocation succeed
       part->n_omega_y= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_y== NULL) printf("n_omega_y: allocation failed"); // check if allocation succeeded
       part->n_omega_x= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_omega_x== NULL) printf("n_omega_x: allocation failed"); // check if allocation succeeded
       part->PA= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->PA== NULL) printf("PA: allocation failed");
       part->DM= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->DM== NULL) printf("DM: allocation failed");
       part->Nb_orb= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Nb_orb== NULL) printf("Nb_orb: allocation failed");
       part->delta= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->delta== NULL) printf("delta: allocation failed");
       part->fomega=(double **)calloc(91,sizeof(double *));
                if (part->fomega== NULL) printf("fomega: allocation failed");
       for(int i=0;i<91;i++) {part->fomega[i]=(double *)calloc(91,sizeof(double)); if (part->fomega[i]==NULL) printf("fomega[i]: allocation failed");}
       part->Smin_fermi= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->Smin_fermi== NULL) printf("Smin_fermi: allocation failed");

       part->temp=(double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->temp== NULL) printf("temp: allocation failed");
       part->Smin_pmps= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin_pmps== NULL) printf("Smin_pmps: allocation failed");
       part->Smin_fast= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Smin_fast== NULL) printf("Smin_fast: allocation failed");
       part->w_int= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_int== NULL) printf("w_int: allocation failed");
       part->w_r_fast= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_r_fast== NULL) printf("w_r_fast: allocation failed"); // check if allocation succeeded
       part->w_r_pmps= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->w_r_pmps== NULL) printf("w_r_pmps: allocation failed"); // check if allocation succeeded
       part->Fx= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fx== NULL) printf("Fx: allocation failed");
       part->sky_XMM= (int *)calloc(part->Npulsars,sizeof(int)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->sky_XMM== NULL) printf("sky_XMM: allocation failed");
       part->sky_chandra= (int *)calloc(part->Npulsars,sizeof(int)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->sky_chandra== NULL) printf("sky_chandra: allocation failed");
       part->n_mu_x= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_mu_x== NULL) printf("n_mu_x: allocation failed");
       part->n_mu_y= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_mu_y== NULL) printf("n_mu_y: allocation failed");
       part->n_mu_z= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_mu_z== NULL) printf("n_mu_z: allocation failed");
       part->ex= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->ex== NULL) printf("ex: allocation failed");
       part->ey= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->ey== NULL) printf("ey: allocation failed");
       part->ez= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->ez== NULL) printf("ez: allocation failed");
       part->nx= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->nx== NULL) printf("nx: allocation failed");
       part->ny= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->ny== NULL) printf("ny: allocation failed");
       part->nz= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->nz== NULL) printf("nz: allocation failed");
       part->cos_in= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->cos_in== NULL) printf("cos_in: allocation failed");
       part->Temp_n= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->Temp_n== NULL) printf("Temp_n: allocation failed");
       part->r_hn= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->r_hn== NULL) printf("r_hn: allocation failed");
       part->PF= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->PF== NULL) printf("PF: allocation failed");
       part->cos_is= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->cos_is== NULL) printf("cos_is: allocation failed");
       part->n_nx= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_nx== NULL) printf("n_nx: allocation failed");
       part->n_ny= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_ny== NULL) printf("n_ny: allocation failed");
       part->n_nz= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->n_nz== NULL) printf("n_nz: allocation failed");
       part->n_sx= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_sx== NULL) printf("n_sx: allocation failed");
       part->n_sy= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->n_sy== NULL) printf("n_sy: allocation failed");
       part->n_sz= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->n_sz== NULL) printf("n_sz: allocation failed");
       part->theta_s= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->theta_s== NULL) printf("theta_s: allocation failed");
       part->theta_n= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->theta_n== NULL) printf("theta_n: allocation failed");
       part->Fxmax= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fxmax== NULL) printf("Fxmax: allocation failed");
       part->Fxmin= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->Fxmin== NULL) printf("Fxmin: allocation failed");
       part->phi_n= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->phi_n== NULL) printf("phi_n: allocation failed");
       part->phi_s= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->phi_s== NULL) printf("phi_s: allocation failed");
       part->mu_hot_ang1= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->mu_hot_ang1== NULL) printf("mu_hot_ang1: allocation failed");
       part->mu_hot_ang2= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->mu_hot_ang2== NULL) printf("mu_hot_ang2: allocation failed");
       part->see_n= (int *)calloc(part->Npulsars,sizeof(int)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->see_n== NULL) printf("see_n: allocation failed");
       part->see_s= (int *)calloc(part->Npulsars,sizeof(int)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->see_s== NULL) printf("see_s: allocation failed");
       part->js= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->js== NULL) printf("js: allocation failed");
       part->jn= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate) 
                if (part->jn== NULL) printf("jn: allocation failed");
       part->r_hs= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->r_hs== NULL) printf("r_hs: allocation failed");
       part->Temp_s= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->Temp_s== NULL) printf("Temp_s: allocation failed");
       part->Lgamma= (double *)calloc(part->Npulsars,sizeof(double)); // (*part->Pinit first elemenet of the table) initialize pointer (allocate)
                if (part->Lgamma== NULL) printf("Lgamma: allocation failed");
}
