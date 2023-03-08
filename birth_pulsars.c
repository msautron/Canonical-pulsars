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




int birth(void *params){//generates Npulsars with an initial period, B and age 

        struct func_params *part= (struct func_params*)params;
	FILE *fp=NULL;
	//char *fname;
        //fname = malloc(50);
	//sprintf(fname,"population_init%.0ld.dat",part->Npulsars);
	/*sprintf(fname,"population_init");
	if ( ( fp = fopen(fname,"w+")) == (FILE *)NULL){
		fprintf(stderr,"Couldn't open file %s \n",fname);
	//exit(-1);
		}*/
        long np=0;

           	while(np<part->Npulsars){  

                        part->Pinit[np]=part->p_mean+gsl_ran_gaussian_ziggurat(part->r,part->sigma_p);
			if (part->Pinit[np] < 0){
		            continue;	
			} 

                        part->Binit[np]=pow(10,log10(part->b_mean)+gsl_ran_gaussian_ziggurat(part->r,part->sigma_b));
	//	        printf("%e %e %e %ld \n",part->Binit[np],part->Pinit[np],part->age_pulsar[np],np);
		        part->age_pulsar[np]  =   part->birth_rate*np*365*24*3600; //s          

/* would be better to print in a file to save time */
 
//		        fprintf(fp,"%e %e %e %ld \n",part->Binit[np],part->Pinit[np],part->age_pulsar[np],np);
		        np++;
                }


//       fclose(fp);
       //gsl_rng_free (part->r);
       return(0);

        }
           
