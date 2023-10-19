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
#include<time.h>

/*

int distribution(void *params){

    	struct func_params *part= (struct func_params*)params;
    	long np=0;
	double two_pi=2*M_PI;   
	double phi,rz;
	double xx,yy,zz; //coordinates relative to eart
//printf ("generator type: %s\n", gsl_rng_name (r));
 // printf ("seed = %lu\n", gsl_rng_default_seed);
 // printf ("first value = %lu\n", gsl_rng_get (r));



           	while(np<part->Npulsars){   

                                                  
	 		phi	       =  two_pi*gsl_rng_uniform(part->r); // galactocentric azimuth uniformly distributed 
                        //rz             =  fabs(gsl_ran_gaussian_ziggurat(part->r,part->sigma_rz)); // gaussian distribution for the distance to the z-axis 
                       // rz             =  gsl_ran_exponential(part->r,part->scale_height); // gaussian distribution for the distance to the z-axis 
                        rz             =  64.6*pow(1.528,4.35)*gsl_sf_gamma(4.35)*gsl_ran_gamma(part->r,4.35,1.528); // distribution taken from Lorimer 2004 
			part->z[np]    =  gsl_ran_exponential(part->r,part->scale_height); //kpc
			if(np%2==0) part->z[np] = part->z[np];
		       	      else part->z[np]  = -part->z[np];	
			part->x[np]   = rz*cos(phi);  //kpc
			part->y[np]   = rz*sin(phi);  //kpc
			xx            = 8.5-part->x[np]; // shift center of the Galaxy to the Sun
			yy            = part->y[np]; //coordinates relative to Earth
			zz	      = part->z[np];	
			part->dist[np]= sqrt(sq(xx)+sq(yy)+sq(zz));
			//printf("phi= %e rGC= %e \n",phi,rGC);
			//printf(" rz= %e x= %e y= %e  z= %e \n",rz,part->x[np],part->y[np],part->z[np]);
			//printf("gsl_ran_gamma(a,b) %e  gamma(4.35) %e pow(a,b) %e \n",gsl_ran_gamma(part->r,4.35,1.258),gsl_sf_gamma(4.35),pow(4.35,1.258));
			//printf(" %e  %e  %e \n",xx,yy,zz);
			//printf("theta %e cos_theta= %e sin_theta= %e cos_phi= %e sin_phi= %e \n",theta,cos_theta,sin_theta,cos(phi),sin(phi));
			np++;
                        //printf("\n");
		}
        //gsl_rng_free (r);


	
	return(0);


   }

*/


int distrib_init(void *params){



        struct func_params *part= (struct func_params*)params;
        long np=0;
        double two_pi=2*M_PI;
	srand((unsigned)time(NULL));
        double phi,R;double p,rng_0_1_d;int rng_0_1;



                while(np<part->Npulsars){
                        phi            =  two_pi*gsl_rng_uniform(part->r); // galactocentric azimuth uniformly distributed 
                        //rz             =  1.0683*pow(4.5,2)*gsl_sf_gamma(2)*gsl_ran_gamma(part->r,2,4.5); // distribution taken from Lorimer 2004 
                        R              =  gsl_ran_gamma(part->r,2,part->Rexp); // distribution taken from packzinski 
                        part->z0[np]    =  gsl_ran_exponential(part->r,part->zexp); //kpc
			rng_0_1=(rand() % 1000000001);rng_0_1_d=rng_0_1;p=rng_0_1_d/1000000000.0;
                        if(p<0.5) part->z0[np] = part->z0[np];
                        else if(p>=0.5)  part->z0[np]  = -part->z0[np];
                        part->x0[np]   = R*cos(phi);  //kpc
                        part->y0[np]   = R*sin(phi);  //kpc
			part->x[np]= part->x0[np];
			part->y[np]= part->y0[np];
			part->z[np]= part->z0[np];
//			printf(" %e %e %e %e \n",R,part->x0[np],part->y0[np],part->z0[np]);
                        np++;


                }

return(0);

}




FILE *kick(void *params){

    	struct func_params *part= (struct func_params*)params;
    	long np;
	double two_pi=2*M_PI;
	double v;
	double cos_theta,phi;
	double vx,vy,vz;
	double x_s,y_s,z_s;
	double dx,dy,dz;
	const double kpc2km=3.0856775807e16; //kpc to km
	//const double kpc2m=3.0856775807e19; //kpc to km
	const double yr_sec=365*24*3600; // year to sec
	double RAD = 180/M_PI;
	double r;
	double glr,gbr; // galactic latitude, galactic longitude in radians
	double age_pulsar_yr;
	FILE *fp_gal=NULL;
	//FILE *file=NULL;
	//int x[100];
	char *fname;
        fname = malloc(150);
	//file=fopen("kick_result.txt","w+");
	if ( ( fp_gal = fopen("galactic_coord.dat","w+")) == NULL){
		fprintf(stderr,"Couldn't open file galactic_coord.dat \n");
	exit(-1);
		}
			
           	for(np=0;np<part->Npulsars;np++){  
  	 		age_pulsar_yr    =  part->age_pulsar[np]/yr_sec; // /!\   should use a local variable
			if (age_pulsar_yr<3*1e6) part->sigma_v=part->v_young;
				else part->sigma_v=part->v_old;

			//Computation of the velocity if we ignore the fact that pulsars are going in the same direction as the rotation axis
			/*vx = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v); //km/s
			vy = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
			vz = gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
        		v  = sqrt(vx*vx + vy*vy + vz*vz); // is  v really useful??*/

			//Computation of the velocity if we do not ignore the fact above
			srand((unsigned)time(NULL));
			cos_theta   =  2*gsl_rng_uniform(part->r)-1;
                        phi         = two_pi*gsl_rng_uniform(part->r);
                        part->n_omega_x[np]=sqrt(1-sq(cos_theta))*cos(phi);
			part->n_omega_y[np]=sqrt(1-sq(cos_theta))*sin(phi);
                        part->n_omega_z[np]=cos_theta;
                        v=sqrt(8.0/M_PI)*part->sigma_v+gsl_ran_gaussian_ziggurat(part->r, part->sigma_v);
                        vx=v*part->n_omega_x[np];
			part->vx[np]=vx;
                        vy=v*part->n_omega_y[np];
			part->vy[np]=vy;
                        vz=v*part->n_omega_z[np];
			part->vz[np]=vz;

			dx = vx*part->age_pulsar[np]/kpc2km;
			dy = vy*part->age_pulsar[np]/kpc2km;
			dz = vz*part->age_pulsar[np]/kpc2km;

                        // coordinates in the Galactocentric system 
			part->x[np]=part->x0[np]+dx; 
			part->y[np]=part->y0[np]+dy;
			part->z[np]=part->z0[np]+dz;
		        	
                        // coordinates relative to the Sun in kpc 
			x_s            = part->x[np];part->x_s[np]=x_s; // shift center of the Galaxy to the Sun
			y_s            = part->y[np]-8.5; part->y_s[np]=y_s;
			z_s            = part->z[np]-0.015; part->z_s[np]=z_s;
    			r	       = sqrt(x_s*x_s+y_s*y_s);

			part->dist[np]= sqrt(sq(x_s)+sq(y_s)+sq(z_s));

			// Galactic coordinates 

			if(part->dist[np]<1e-15){
			     gbr=0;
			 }else{
			     gbr=asin(z_s/part->dist[np]);
			 }
			part->gb[np]=gbr*RAD;
			    
			if(r<1e-15){
			      glr=0;
			}else{
			      if(x_s>=0){
				glr=acos(-y_s/r);
			      }else{
				glr=acos(y_s/r)+M_PI;
			      }
			 }

			part->gl[np]=glr*RAD;

			//part->gb[np] = RAD*asin(part->z[np]/part->dist[np]);
			//part->gl[np] = 180*acos(xx/part->dist[np]/M_PI);
			//part->gb[np] = 180*asin(part->z[np]/part->dist[np]/M_PI);
			//part->gb[np] = 180*asin(part->z[np]/part->dist[np])/M_PI;


			fprintf(fp_gal,"Gal %e %e %e %d \n",part->gl[np],part->gb[np],part->dist[np],2);
			//fprintf(file,"%e %e\n",part->x[np],part->y[np]);
	//		printf("Gal %e %e %e %d \n",part->gl[np],part->gb[np],part->dist[np],2);
			//printf("%e \n",part->gl[np]);
			//printf("%e  %e \n",part->x[np],part->z[np]);
        	} 
		//fclose(file);


return fp_gal;

/*
    FILE *fp=NULL;
       if ( ( fp= fopen("file.dat","w+")) == NULL){
                printf("Couldn't open file file.dat \n");
                exit(-1);
        }
       int i;
    for (i=0;i<10;i++){
	x[i]=i*i;	
 	//printf("%d \n", x[i]);
 	fprintf(fp,"%d \n",x[i]);
 //	fprintf(fp,"%s \n", "This is file.dat");
    }
    return fp;
*/
}



int geometry(void *params){ // calculated the geometry of the pulsar (i.e xi, and rho) 


    	struct func_params *part= (struct func_params*)params;
    	long np=0;
	//double two_pi=2*M_PI;   
	//double phi,cos_theta;
	double Npulsars=part->Npulsars;
	double n[3];
	//double n_omega[3];
	double norm;
	double xi,alpha;
	double rho;
	double ratio;
	double k=0.118940963;
	
	
           	for (np=0;np<Npulsars;np++){ 

			part->np       = np;
                        alpha          = part->alpha[np]; // just to make it easier to read                         
			//rho            = 3*sqrt(M_PI*hem/(2*part->period[np]*SI_C)); // replace by a numeric value
			rho            = k/sqrt(part->period[np]); // 3*sqrt((PI*hem)/2*P*c) Equation 2 of Johnston et al. (2020)
			part->rho[np]  = rho;
	 		//phi	       =   two_pi*gsl_rng_uniform(part->r); // azimuth of the rotation axis 
			//if(np%2==0)  cos_theta   =  gsl_rng_uniform(part->r); // theta of the rotation axis
		       	  //    else   cos_theta   =  -gsl_rng_uniform(part->r); // theta is between -1 and 1 
			//cos_theta   =  gsl_rng_uniform(part->r); // theta of the rotation axis

			/* line of sight vector */
			norm=sqrt(sq(part->x[np])+sq(part->y[np]-8.5)+sq(part->z[np]-0.015));
			n[0]=(part->x[np])/norm;  
			n[1]=(part->y[np]-8.5)/norm;
			n[2]=(part->z[np]-0.015)/norm;

			/* unit vector of the rotation axis */
			/*n_omega[0]=((1-sq(cos_theta))*cos(phi));part->n_omega_x[np]=n_omega[0];
			n_omega[1]=((1-sq(cos_theta))*sin(phi));part->n_omega_y[np]=n_omega[1];
			n_omega[2]=cos_theta;part->n_omega_z[np]=n_omega[2];*/

			/* angle between vectors n and n_omega	*/
		        xi=acos(part->n_omega_x[np]*n[0]+part->n_omega_y[np]*n[1]+part->n_omega_z[np]*n[2]); 
			part->xi[np]=xi;



			/* calculates the width of the radio beam w_r, from eq 22 of our paper */ 
				if(fabs(alpha-xi)<= rho && alpha >= rho && xi < M_PI/2){ 
					ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
					part->w_r[np]=2*acos(ratio);
				} else if(fabs(xi-(M_PI-alpha))<=rho && alpha >= rho && xi > M_PI/2){ 
					alpha=M_PI-alpha;
					ratio=(cos(rho)-cos(alpha)*cos(xi))/(sin(alpha)*sin(xi));
					part->w_r[np]=2*acos(ratio);		
				}	
		}

    return(0);
}


void spinvel_angle(void *params){

	struct func_params *part= (struct func_params*)params;
	long np=0;
	double norm_v;double norm_omega;
	for(np=0;np<part->Npulsars;np++){
	   norm_v=sqrt(sq(part->vx[np]*1e3)+sq(part->vy[np]*1e3)+sq(part->vz[np]*1e3));
	   norm_omega=sqrt(sq(part->n_omega_x[np])+sq(part->n_omega_y[np])+sq(part->n_omega_z[np]));
	   part->PA[np]=(part->n_omega_x[np]*(part->vx[np]*1e3)+part->n_omega_y[np]*(part->vy[np]*1e3)+part->n_omega_z[np]*(part->vz[np]*1e3))/(norm_v*norm_omega);
           //if (np==0) printf("value computed for PA: %e",part->PA[np]);
	}

}

int detection(void *params){ //check the flux of each pulsar and if the beam sweps the Earth


    	struct func_params *part= (struct func_params*)params;
    	long np=0;
	double B;double P;
	long count_radio_tot=0;
	long count_radio=0;
	long count_gamma=0;
	long count_radio_gamma=0;
	long Nsup_gamma=0;
	long Nsup2_gamma=0;
	long Nsup=0;
	long Nsup2=0;
	long Nsup_radio=0;
	long Nsup_radio_gamma=0;
	long Nsup2_radio=0;
	long Nsup2_radio_gamma=0;
	long Nr=0;
	long Ng=0;
	double xi,rho,alpha;
	double Smin_gamma, Smin_radio;
        double glat;
	//double ratio;
	long Ntot;
	long Nbeam=0;
        int detec;
	double S_N ; // signal to noise ratio
	long S_Nmin=10 ; // signal to noise ratio
	//FILE *fp=NULL;
	FILE *file_data=NULL;
	FILE *check_val;
	//FILE *check_val2;
	double eta=0.15;double alpha_l=45*(M_PI/180);double T6=2;double b=40;
	double alpha_l2;double T6_2;
	double P_dot_line;
	//char *fname;
        //fname = malloc(150);
	file_data=fopen("P_Pdot_positions.txt","w+");
	check_val=fopen("check_val.txt","w+");
	//check_val2=fopen("check_val2.txt","w+");
			
           	for (np=0;np<part->Npulsars;np++){ 
	               // np=part->np;
			xi=part->xi[np];
			B=part->B[np];
			P=part->period[np];
			//ratio=B/sq(P);
			rho=part->rho[np];
			alpha=part->alpha[np];
			P_dot_line=(3.16e-4*pow(T6,4)*sq(part->period[np])*1e-15)/(sq(eta)*b*sq(cos(alpha_l)));
			if (P_dot_line/part->Pdot[np] >= pow(10,-0.4) && P_dot_line/part->Pdot[np] <= pow(10,0.4)) {
				alpha_l2=20*(M_PI/180)*gsl_rng_uniform(part->r);
				T6_2=6*gsl_rng_uniform(part->r)+1;
				P_dot_line=(3.16e-4*pow(T6_2,4)*sq(part->period[np])*1e-15)/(sq(eta)*b*sq(cos(alpha_l2)));
			}
			glat = (180*asin(part->z[np]/part->dist[np]))/M_PI;
			//glat = 180*asin(part->z[np]/part->dist[np]/M_PI);
			//wtilde=part->w_r[np]*part->period[np]/(2*M_PI); // converting the width from radians to s 
			//Smin_gamma = 5e-12;
			Smin_radio = part->Smin[np];
		//	Smin_radio = 0.15;
			//S_N = part->Fr[np]/part->Smin[np];
			//S_N = part->Fr[np]/Smin_radio;
			if (part->NenuFAR==1 ) S_N = part->flux_low_freq[np]/Smin_radio; // nenuphar
			else S_N = part->Fr[np]/part->Smin[np];

                        if (S_N > S_Nmin) detec = 1;
				else detec  = 0;
                 


			if ((((-21.80 < glat && glat < 27.13) && part->NenuFAR==1)) || part->NenuFAR==0){

			        					        
			/* radio detection */
				if(fabs(alpha-xi)<= rho && alpha >= rho){ 
					Nbeam++;
                			if (detec==1){
						Nr=1;  // mJy
						part->detec[np]=1;
                                                part->detec_rad[np]=1;
						count_radio_tot++;
						//fprintf(check_val2,"%e|%e|%e\n",B,P,P_dot_line);
					}
				} else if (fabs(xi-(M_PI-alpha))<=rho && alpha >= rho){ 
					Nbeam++;
                			if (detec==1){
						count_radio_tot++;
						part->detec[np]=1;
                                                part->detec_rad[np]=1;
						Nr=1;  
						//fprintf(check_val2,"%e|%e|%e\n",B,P,P_dot_line);
						}
					}
				}
			
			/* gamma detection */
			  	if(fabs(xi-M_PI/2.)<=alpha){ 
					if(Nr==1 && glat < 2.) {  //deep follow-up surveys at 1.4 GHz of gamma-ray sources
					    Smin_gamma=4e-15; // in W.m^-2 
					} else Smin_gamma=16e-15;
				//Smin_gamma=5e-12;
                			if(part->Fg[np]>Smin_gamma) {Ng=1;part->detec[np]=1;part->detec_gam[np]=1; }  
				}  else continue; 


				if(Ng==1 && Nr==1){  //both radio and gamma are detected
				        if (P_dot_line>part->Pdot[np]) {Ng=0;Nr=0;part->detec_rg[np]=0;part->detec[np]=0;part->detec_rad[np]=0;part->detec_gam[np]=0;}
					else{
					   count_radio_gamma++;part->detec_rg[np]=1;part->detec[np]=1;part->detec_rad[np]=0;part->detec_gam[np]=0;	
		
			        	       	   if(part->Edot[np]>1e28) { //28
							   Nsup_radio_gamma+=1;
						   }
			        		   if(part->Edot[np]>1e31) { //31
							   Nsup2_radio_gamma+=1;
						   }
						   fprintf(check_val,"%e|%e|%e\n",B,P,P_dot_line);
				} 
				}
				else if(Nr==1 && Ng==0){ //radio only
                       
					if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_rad[np]=0;}
					else{
					   count_radio++;part->detec[np]=1;part->detec_rad[np]=1;

			        		   if(part->Edot[np]>1e28) { 
							   Nsup_radio+=1;
						   }
			        		   if(part->Edot[np]>1e31) { 
							   Nsup2_radio+=1;
						   }
						   fprintf(check_val,"%e|%e|%e\n",B,P,P_dot_line);
				}
				}
				else if(Ng==1 && Nr==0){ //gamma only
				
					if (P_dot_line>part->Pdot[np]) {part->detec[np]=0;part->detec_gam[np]=0;}
					else{
					   count_gamma++;part->detec[np]=1;part->detec_gam[np]=1;

			        		   if(part->Edot[np]>1e28) { 
							   Nsup_gamma+=1;
						   } 
			        		   if(part->Edot[np]>1e31) { 	
							   Nsup2_gamma+=1;
						   }
						   fprintf(check_val,"%e|%e|%e\n",B,P,P_dot_line);
				} 
				}

				if(Ng==1 || Nr==1) {
			        		if(part->Edot[np]>1e28) { 
							Nsup++;
						}
			        		if(part->Edot[np]>1e31) { 
							Nsup2++;
						}	
				}
			Nr=0;
			Ng=0;
			detec=0;
			       if(part->detec_rad[np]==1){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|1|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->cos_a0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np]);
                        }

                               else if(part->detec_gam[np]==1){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|2|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->cos_a0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np]);
                        }

                               else if(part->detec_rg[np]==1){

                                  fprintf(file_data,"%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|3|\n",part->period[np],part->Pdot[np],part->x[np],part->y[np],part->age_pulsar[np],part->err_rel_g[np],part->dist[np],part->gl[np],part->gb[np],part->cos_a0[np],part->alpha[np],part->B[np],part->z[np],part->vx[np],part->vy[np],part->vz[np],part->vx0[np],part->vy0[np],part->vz0[np],part->PA[np]);
                        }


		}

 
			Ntot = count_radio+count_gamma+count_radio_gamma;


			printf("# Number of pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup);
			printf("# Number of radio-only pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_radio);
			printf("# Number of gamma-only pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_gamma);
			printf("# Number of radio+gamma pulsars with Edot > 1e32 erg.s-1: %ld \n",Nsup_radio_gamma);
			printf(" \n");
			printf("# Number of pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2);
			printf("# Number of radio-only pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_radio);
			printf("# Number of gamma-only pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_gamma);
			printf("# Number of radio+gamma pulsars with Edot > 1e33 erg.s-1: %ld \n",Nsup2_radio_gamma);
			printf(" \n");
			printf("# Total number of detected pulsars  %ld \n",Ntot);
			printf("# Total number of radio-only pulsars detected %ld \n",count_radio);
			printf("# Total number of gamma-only pulsars detected %ld \n",count_gamma);
			printf("# Total number of radio+gamma pulsars detected %ld \n",count_radio_gamma);
			printf("# Number of radio pulsars beaming to us %ld \n",Nbeam);
			printf(" \n");
			fclose(file_data);
			fclose(check_val);
			//fclose(check_val2);

return(0);
}

/* /!\  units should be changed to SI units */
int radio_flux(void *params){ // returns Smin and radio flux Fr in mJy 
    	struct func_params *part= (struct func_params*)params;
   	double Fj; //scatter term
	long np=0;
	double wtilde;
	double S_0=0.05;

           	for (np=0;np<part->Npulsars;np++){ 
		
			part->np=np;	
			wtilde=part->w_r[np]*part->period[np]/(2*M_PI); // converting the width from radians to s 


		//	 	part->Smin[np]=8.660254e-05*sqrt(wtilde/((part->period[np])-wtilde)); // Bhattcharya 1998 this if for S/N=10, BW=100 MHz, G=10K/Jy, tau=10 min, Tsys= 30K, np=2  
		  //  		part->Smin[np]=part->Smin[np]*1e3; // converts from Jy to mJy
		//	} else {
		        if(part->NenuFAR==1) S_0=30.;

			else if (part->ska==1) S_0 = 0.0045;

			else S_0 = 0.05;

		        part->Smin[np]=S_0*sqrt(wtilde/((part->period[np])-wtilde)); //in mJy	
			//printf("Smin %e \n",part->Smin[np]);
				 //part->Smin[np]=15e-3; 	

                        
			//part->Smin[np]=0.15;
			//part->Smin[np]=0.1;
			Fj             =  fabs(gsl_ran_gaussian_ziggurat(part->r,0.2)); 
                    	part->Fr[np]   =(9.0/sq(part->dist[np]))*(pow((part->Edot[np]*1e7),0.25)/1e9)*pow(10,Fj); // (From Johnston's paper) in mJy .1e7 to get Edot in erg.s-1
			//printf("Fr %e dist %e Edot %e \n",part->Fr[np],part->dist[np],part->Edot[np]);
			//printf("S_0 %e Smin %e Fr %e sqrt %e wtilde %e P %e \n",S_0,part->Smin[np],part->Fr[np],sqrt(wtilde/((part->period[np])-wtilde)),wtilde,part->period[np]);
//			printf("S_0 %e Smin %e sqrt %e wtilde %e P %e \n",S_0,part->Smin[np],sqrt(wtilde/((part->period[np])-wtilde)),wtilde,part->period[np]);
			// /!\ problem: units !!!
			//printf("Fr %e np %ld dist %e \n",part->Fr[np],np,part->dist[np]);		
			//printf("%e \n",part->Fr[np]);		

		}

 



	return(0);

}



int gamma_flux(void *params){

    	struct func_params *part= (struct func_params*)params;
   	double cg,fac; //scatter term
	cg=1.;
	double lgamma;
	const double kpc2m=3.0856775807e19; //kpc to m
	fac=1./(4*M_PI*cg);
	//fac=0.3;
	long np=0;

           	for (np=0;np<part->Npulsars;np++){ 
		
			part->np=np;	

                        //lgamma = 1e15*pow(part->B[np]*1e4,0.11)*pow(part->Edot[np]*1e7,0.51); 
                        lgamma = pow(10,26.15)*pow(part->B[np]/1e8,0.11)*pow(part->Edot[np]/1e26,0.51); //(W) from Kalapotharakos et al 2019
		    
	        	if(part->alpha[np]<-part->xi[np]+0.6109) cg=1.9; //the correction factor cg (or f_omega)
		   		else cg=1.;
      			part->Fg[np]=fac/(sq(part->dist[np]*kpc2m))*lgamma; // gamma flux in W.m-2
			//printf("Fg %e np %ld dist %e \n",part->Fg[np],np,part->dist[np]);		
                        //printf("lgamma %e Fg %e B %e Edot %e \n",lgamma,part->Fg[np],part->B[np],part->Edot[np]);	

      			//part->Fg[np]=fac*(part->Edot[np]*1e7/(sq(part->dist[np]*KPC2CM)))*sqrt(1e33/(part->Edot[np]*1e7)); // Petri 2011a 
 				
		} 

return(0);


}

/*void verif_log_normal(void *params){

	struct func_params *part= (struct func_params*)params;
	FILE *file=NULL;
	file=fopen("verif_lognorm.txt","w+");
	for(int i=0;i<part->Npulsars;i++){
           fprintf(file,"%e\n",part->Binit[i]);
	}
	fclose(file);

}*/

double calc_Pdotline(void *params,long np){

        struct func_params *part= (struct func_params*)params;
        double log_pdot_line;
        log_pdot_line=3*log10(part->period[np])+log10((16*(pow(M_PI,3))*pow(R_NS,6)*1.5)/(SI_I*SI_mu0*pow(SI_C,3)))+2*(log10(0.17)+8.0);
        return log_pdot_line;
}

double is_dead(double B,double P){ // Death line from Faucher Giguère & Kaspi (2006) 

	//struct func_params *part= (struct func_params*)params;
	double ratio=B/sq(P);
	return ratio;
	//double ratio2=part->Pdot[np]/(part->period[np]*sq(part->period[np]));
	//double cmp_value=3*log10(part->period[np])+log10(16*sq(M_PI)*M_PI*sq(R_NS)*sq(R_NS)*sq(R_NS)*(1.0+sq(sin(45*M_PI/180)))*sq(0.17e8)/(SI_I*SI_mu0*sq(SI_C)*SI_C));
	//if (log10(part->Pdot[np]) >= cmp_value) {return 1;}
	//else if (log10(part->Pdot[np]) < cmp_value) {return 0;}
	//if (ratio < 0.17e8) {return 0;}
	//else if (ratio > 0.17e8) {return 1;}
}

/*int is_dead2(void *params,long np){ // Death line from Beskin & Istomin 2022 -> too restrictive 

	struct func_params *part= (struct func_params*)params;
	double xsi=11;double lambda=41;double f_star=1.9;double Kg=0.07;double F=0.7;double h_x0=3.1;
	double beta_d=2.1*pow(xsi,-1/2)*Kg*pow((f_star/1.6),-9/4)*pow(lambda/15,-1/2)*h_x0*F;
	double P_dot_death=1e-15*beta_d*pow(part->period[np],11/4);
	if (part->Pdot[np]<P_dot_death) return 0;
	else if (part->Pdot[np]>=P_dot_death) return 1;

}*/

/*int is_dead3(void *params,long np){ // Death line from Wu et al. (2020) -> results similar to death line one 

        struct func_params *part= (struct func_params*)params;
        if (part->Edot[np] < 1e24) {return 0;}
        else if (part->Edot[np] > 1e24) {return 1;}
}*/

int is_dead4(void *params,long np){ // Death line from Mitra et al. (2019) 

        struct func_params *part= (struct func_params*)params;
	double eta=0.15;double alpha_l=45*(M_PI/180);double T6=2;double b=40;
	double P_dot_lim=(3.16e-4*pow(T6,4)*sq(part->period[np])*1e-15)/(sq(eta)*b*sq(cos(alpha_l)));
	if (part->Pdot[np]>=P_dot_lim) return 1;
	else return 0;
}

int radio_flux_low_freq(void *params){

    	struct func_params *part= (struct func_params*)params;
	long np=0;
	double S_0;
	double alpha_positiv;
	double alpha_negativ;
	double nu_break;
        double scale_pos=0.812;
	double shape_pos=0.873;
	double loc_pos=-0.014;
	double shape_break=0.698;
        double scale_break=152.9;
	double loc_break=0.0;
	double shape_neg=0.2;
        double scale_neg=3;
	double loc_neg=-4.6;
	double S_1PL,S_2PL;
	double nu_0=1.4e9; //1.4 GHz
	double nu_nenuphar = 58e6 ; 
	double nu=nu_nenuphar;
 
	//nu= 1e7;	
	np=0;
           	for (np=0;np<part->Npulsars;np++){ 
        //   	while(np<part->Npulsars){  
			part->np=np;	
	//		np = part->np;
 			S_0           =  part->Fr[np];
		//	wtilde        =  part->w_r[np]*part->period[np]/(2*M_PI); // converting the width from radians to s 
                    	alpha_positiv=  gsl_ran_lognormal(part->r,log(scale_pos),shape_pos); 
                        alpha_positiv =  (alpha_positiv-loc_pos);
			//printf("wtilde %e period %e \n",wtilde,part->period[np]);

                    	nu_break      =  gsl_ran_lognormal(part->r,log(scale_break),shape_break); // in MHz
			nu_break      =  (nu_break-loc_break);
			nu_break      =  nu_break*1e6;
                    	alpha_negativ =  gsl_ran_lognormal(part->r,log(scale_neg),shape_neg); // Problem: Fj is the same each time the function radio_flux is called...
                        alpha_negativ =  (alpha_negativ+loc_neg);
           	//	while(nu<5*1e10){  
			/*if (nu_break > 1.4e9){
		            continue;	
			}
*/
					 //printf("%e %e \n",nu,S_2PL);
				//	 printf("# nu_break %e alpha_pos %e alpha_neg %e \n",nu_break,alpha_positiv,alpha_negativ);
 //         alpha_negativ = -1.7;
//	S_0=100;
			//printf("%e %e \n",nu,part->flux_low_freq[np]);
			//nu = pow(10,lognu);
		//	printf("%e %e \n",nu,lognu);
         //  	        while (nu<5e10){ 


				if (np%2 == 0){

					 S_1PL                   = S_0*pow((nu/nu_0),alpha_negativ);
					 part->flux_low_freq[np] = S_1PL;	
				//	 printf("%e %e %e \n",nu,S_1PL,part->Fr[np]);
					 //printf("%e %e \n",nu,S_2PL);
					 //printf("Flux (1PL), nu < nu_break %e S_1400 %e nu_break %e alpha_pos %e alpha_neg %e \n",S_1PL,S_0,nu_break,alpha_positiv,alpha_negativ);
                                } 
				else {
					if(nu<nu_break) {

					 S_2PL                 = S_0*pow((nu_0/nu_break),alpha_positiv)*pow((nu_break/nu_0),alpha_negativ)*pow(nu/nu_0,alpha_positiv);
					 part->flux_low_freq[np] = S_2PL;
				//	 printf("%e %e %e \n",nu,S_2PL,part->Fr[np]);
					 //printf("Flux (2PL), nu < nu_break %e S_1400 %e nu_break %e alpha_pos %e alpha_neg %e \n",S_2PL,S_0,nu_break,alpha_positiv,alpha_negativ);

					}

					else {

					 S_2PL                   = S_0*pow((nu/nu_0),alpha_negativ);
					 part->flux_low_freq[np] = S_2PL;				
					 //printf("%e %e \n",nu,S_2PL);
				//	 printf("%e %e %e \n",nu,S_2PL,part->Fr[np]);
					 //printf("Flux (2PL), nu < nu_break %e S_1400 %e nu_break %e alpha_pos %e alpha_neg %e \n",S_2PL,S_0,nu_break,alpha_positiv,alpha_negativ);
                                        }
			
				}
                            
			//printf("%e \n",nu_break);
			//printf("%e \n",alpha_negativ);
			//printf("%e \n",alpha_positiv);
	       //               }
                       
		//	printf("nu = %e Flux %e alpha_neg %e Flux (1.4 GHZ) %e \n",nu,part->flux_low_freq[np],alpha_negativ,S_0);
			//printf("nu = %e Flux %e alpha_neg %e nu_break  %e \n",nu,part->flux_low_freq[np],alpha_negativ,nu_break);
	//		printf("nu = %e Flux %e alpha_pos %e alpha_neg %e nu_break  %e S_0 %e \n",nu,part->flux_low_freq[np],alpha_positiv,alpha_negativ,nu_break,S_0);
	//		printf("log10(nu/nu_0) %e \n",log10(nu/nu_0));
	//		printf("nu %e Flux %e Smin %e \n",nu,part->flux_low_freq[np],part->Smin[np]);
// 			nu *= 1.1; 
                     np++;

               } 

return(0);
}
