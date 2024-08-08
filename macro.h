#ifndef _MACRO_DEFINED
#define _MACRO_DEFINED
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<stdbool.h>
#define cube(a) ((a)*(a)*(a))
#define maximum(a,b) ((a)>(b)?(a):(b))
#define minimum(a,b) ((a)<(b)?(a):(b))
#define N_MAX   10000000 // maximum number of pulsars
#define sq(a) ((a)*(a))

struct func_params{

	double *n_omega_x;
	double *n_omega_y;
	double *n_omega_z;
	double *PA;
        long birth_rate;
        double R; //Radius of the neutron star
        double *alpha;//angle between the magnetic field and the rotation axis
        double *Binit; //table which returns the initial B
        double *B; //table which returns the initial B
        double *Pinit; //same for P
	double p_mean; // the standard deviation for the initial distribution of p
        double sigma_p; //standard deviation for the distribution of the period
        double b_mean; // idem for b
        double sigma_b;//idem
        double *age_pulsar; // stores the age of the pulsar (in s)
	double *DM; //Stores the dispersion measure of the pulsar in cm^-3.pc
	double *Nb_orb; //Estimation of the number of orbits done by the NS 
	double *S_N;
	double Bfield_const;
	long Bfield_var;
        double *x; //kpc  x coordinate in the Galctocentric frame
        double *y;//idem 
        double *z; // kpc height to the galactic plane
	double *gb; // galactic latitude
	double *gl; // galactic longitude
        double *z0; //kpc initial z0 in the Galactocentric frame
        double *x0; // idem 
        double *y0; //idem
        double *x_s;
	double *y_s;
	double *z_s;
	double *vx0; //table which stores the initial velocities on the x absciss, km/s
	double *vy0; //table which stores the initial velocities on the y absciss, km/s
	double *vz0; //table which stores the initial velocities on the z absciss, km/s
        double *vx;  // velocities on the x absciss, km/s
	double *vy;  // velocities on the y absciss, km/s
	double *vz;  // velocities on the z absciss, km/s
        double *err_rel_g; // relative error on the energy for the integration of the equation of movement 
        long   np; //pulsar number= index of tables that store the pulsars parameters) 
        double sigma_v; // sigma 1D
        double *period; //stores the actual period
        long    Npulsars; //total number of pulsar
	double Rexp; //in kpc parameter for the R distribution
	long v_young; //sigma 1D for the young pulsarvelocity
	long v_old; // old
        long *detec; // indicates if the pulsar is considered as detected
	long *detec_rad; // indicates if the pulsar is detected in radio
	long *detec_gam; //indicates if the pulsar is detected in gamma
	long *detec_rg; // indicates if the pulsar is detected in radio and gamma
        double zexp; // in kpc
	double *Smin;
	double Fmin;
	double Kr;
	double *Edot;
	double *Pdot;
	double *xi; //angle in radians
	double *rho; //witdh of the beam in radian
	double *w_r; //width of the radio profile
	double *Fr; //radio flux table
	double *flux_low_freq; //radio flux table
	double *Fg; //gamma
        double *cos_a0; //a0= initial inclination angle
 	double tau_MHD_al;
 	double tau0_B0;
	double tau0_B0_2;
	double tau0_B0_3;
	//double tau0_B0_4;
	long k_tau0_B0;
 	double alpha_d; // dB/dt = -a B^(1+alpha_d)
	double tau_vac_al;
	double tau_d;
	long ska;
	double *dist; //disance from us
        _Bool vacuum,vacuum_evol,ff_evol;	
        _Bool NenuFAR;	
	const gsl_rng_type * T;
    	gsl_rng * r;  /* global generator */
};

#define SI_C 2.997924858e8 /* speed of light in units of m/s */
#define TMILKY 13.5e9 /* Age of the Milky Way in years */
#define KPC2CM 3.0856775807e21 /* kpc in cm */
#define SI_I 1e38   /*Moment of Inertia in kg.m2 */ 
#define R_NS 12000   /*Moment of Inertia in kg.m2 */ 
#define M_PI 3.14159265358979323846 /* pi */
#define SI_mu0 1.25663706212e-6 /* vacuum permeability in H/m */
#define G_grav 6.67430e-20 //gravitational constant km^-3 kg^-1 s^-2
#define MSUN 1.98847e30 //Solar mass kg
#endif

