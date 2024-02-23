#include"macro.h"
#include <stdlib.h>
#include <stdio.h>
#include"initialize.h"
#include<gsl/gsl_rng.h>
#include <math.h>
#include "birth_pulsars.h"
#include "detection.h"
#include "evolution.h"
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf_gamma.h>
#include<string.h>

#define TAILLE_MAX 1000 // Tableau de taille 1000



int get_dispersion_measure(void *params){

        struct func_params *part= (struct func_params*)params;

 	  /* use the YMW16 model to get the DM measure as a function of gl, gb and d */
	int Npulsars = part->Npulsars;
	FILE *fp=NULL; //stores the output of popen command (./ymw16 ...) in a file
	char path_output[1035];//stores the output of popen command (./ymw16 ...) in a string
	char path_dmtau[1035];//stores the output of popen command (./ymw16 ...) in a string
	FILE *in_file=NULL; // pointer to the file galactic_coord.dat
	FILE *output_file;
	FILE *output_file_dmtau;
	FILE *Extract_dmtau=NULL;
	char fname[1000]="StringInit"; // stores the content of galactic_coord.dat
	in_file = kick(part); //file that contains gl gb dist 2


//	if ( ( in_file = fopen("galactic_coord.dat","r")) == NULL){ //why casting the NULL pointer??????
	if (in_file == NULL){ //why casting the NULL pointer??????
		printf("Couldn't open file galactic_coord.dat \n");
          	exit(-1);
	}
        fseek(in_file, 0, SEEK_SET);

	if ( ( output_file = fopen("output_dispersion_measure.dat","w+")) == NULL){
		printf("Couldn't open file output_dispersion_measure.dat \n");
		exit(-1);
	}

	 //  reads galactic_coord.dat line by line and stores stream in fname 	

          //printf("coucou"); 
//		fp = fopen("test_fscanf.txt","r");

//	fgets(fname,1000,in_file);
	//printf("adresse de in_file %p \n",in_file);
	//printf(" %p %p %p \n",in_file,fgets(fname,1000,in_file),in_file);
	//printf(" fname = %s \n",fname);
	 while (fgets(fname,1000, in_file) != NULL) { //reads data until it either encounters a newline of exceeds specified buffer length
		 printf("fname %s ",fname);
         
  //             char command[1000]="ls -l";
                 char command[1000]="./ymw16 ";
		 strcat(command,fname);  // concaténation of ./ymw16 and Gal gl gb dist 2 for each line of galactic_coord.dat file
		 //printf("après la concaténation %s \n",command);
		 // Open the command for reading. 
		 fp = popen(command, "r");
		 if (fp == NULL) {
		   printf("Failed to run command\n" );
		   exit(1);

		 }

		  // Reads the output a line at a time - output it. 
		 while (fgets(path_output, sizeof(path_output), fp) != NULL) {
		//   	 printf("output is: %s \n", path_output);
		         printf("path_output %s \n",fgets(path_output, sizeof(path_output), fp)); // why is it printing NULL pointer???????
			 fprintf(output_file,"%s ", path_output);
		 } 


		 pclose(fp);

	 }

	 fclose(output_file);
	 
         Extract_dmtau = popen("gawk '{print $9, $11}' output_dispersion_measure.dat", "r");

	 if (Extract_dmtau == NULL) {
	   printf("Failed to run command\n" );
	   exit(1);
	 }

        FILE *pixel = fopen("pixel.txt", "r+");
	if (pixel == NULL)   // check if file could be opened
	{
	printf("Can't open file");
	exit(1);
	}
  	double dm[Npulsars], tau[Npulsars]; 
  	//int arr[Npulsars], arr2[Npulsars]; 

	if ( ( output_file_dmtau = fopen("dmtau.dat","w+")) == NULL){
		printf("Couldn't open file dmtau.dat \n");
		exit(-1);
	} 
/*	if ( ( output_file_dmtau = fopen("dmtau_bckp.dat","r+")) == NULL){
		printf("Couldn't open file dmtau.dat \n");
		exit(-1);
	}*/
	 // Now reads data from dmtau.dat and stores it in arrays dm[], tau[]  
	/*if (dmtau == NULL)   // check if file could be opened
	{
            printf("Couldn't open file dmtau.dat \n");
	    exit(1);
	}
*/
	 while (fgets(path_dmtau, sizeof(path_dmtau), Extract_dmtau) != NULL) {
//	 fgets(path_dmtau, sizeof(path_dmtau), Extract_dmtau) != NULL;  
		 printf("output is dmtau: %s ", path_dmtau);
		 fprintf(output_file_dmtau,"%s ",path_dmtau);
	 } 

	pclose(Extract_dmtau);
	int nbofvaluesread = 0;
	int readitem = 0;

		//for(int i = 0; i < Npulsars; i++)  // read 40000 values
		for(int i = 0; i < 4 ; i++)  // read 40000 values
		{
			//readitem = fscanf(output_file_dmtau,"%le %le",&dm[i],&tau[i]);
			//readitem = fscanf(pixel,"%d %d",&arr[i],&arr2[i]);
   			readitem = fscanf(output_file_dmtau,"%lf %lf",&dm[i],&tau[i]);
   			printf("readitem %d \n",readitem);
   			 //printf("Npulsars %d \n",Npulsars);
		         //printf("dm %e tau %e \n", dm[i],tau[i]);
	//		printf("i =  %d \n",i);
		//	if (fscanf(dmtau,"%lf %lf", &dm[i],&tau[i]) != 2){
			//readitem = fscanf(output_file_dmtau,"%lf %lf", &dm[i],&tau[i]);
                        if (readitem!=2) {
			     printf("could not read file \n");	     
   			     printf("readitem %d \n",readitem);
			     break;   // stop loop if nothing could be read or because there
			}	  					// are less than Npulsars values in the file, or some 

   			     printf("readitem %d \n",readitem);
		  //           printf("dm %lf tau %lf \n", dm[i],tau[i]);
									
			nbofvaluesread++;
	                printf("i %d dm[i] %lf tau[i] %lf number of values read %d \n",i,dm[i],tau[i],nbofvaluesread);
		} 

		printf(" nbof %d \n",nbofvaluesread);
/*
		for(int i = 0; i < nbofvaluesread ; i++){
		   printf("dm %e tau %e \n", dm[i],tau[i]);
		}	
*/

	 /* close */
//	 pclose(Extract_dmtau);
	fclose(output_file_dmtau);
//	 pclose(dmtau);
	 //pclose(fp2);
	 fclose(in_file);
	 //pclose(fp);


return(0);

}

/* result of the convolution of the intrinsic pulse (Gaussian shape) */
/*
double convolution(void *params){

            
    	 struct func_params *part= (struct func_params*)params;
	 double convol; 
	 double sigma;
	 double tau; // dispersion measure
	 long np;
	 double S; // integrated flux density of the pulse at 58 MHz
	 double t;

                for (np=0;np<part->Npulsars;np++){ 
	 		S      = part->flux_low_freq[np]; // integrated flux density of the pulse at 58 MHz
	 		convol = (S/2./tau)*exp(sq(sigma)/2./sq(tau))*exp(-t/tau)*(1+erf((t-sq(sigma)/tau)/sigma*sqrt(2)));



		}

return(convol);


}


*/
