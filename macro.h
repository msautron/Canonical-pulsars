#ifndef _MACRO_DEFINED
#define _MACRO_DEFINED
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<stdbool.h>
#define cube(a) ((a)*(a)*(a))
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define N_MAX   10000000 // maximum number of pulsars
#define sq(a) ((a)*(a))

struct func_params{

	double M_for_K; //Mass of neutron stars used for the computation of compactness
	double R_for_K; //Radius of neutron stars used for the computation of compactness
	double A_propto; //Factor for relation between T and P, Pdot <<north>>
	//double A_propto2; //Factor for relation between T and P, Pdot <<south>>
	double D_propto; //Factor for relation between size of hotspot and sqrt(R_NS/R_L) <<north>>
	//double D_propto2; //Factor for relation between size of hotspot and sqrt(R_NS/R_L) <<south>>
	//double pcst; //Power for the gamma luminosity law (constant)
	//double pb; //Power for the gamma luminosity law (Magnetic field)
	//double pe; //Power for the gamma luminosity law (Spin down luminosity)
	long birth_rate; //BR for the number thres1 of pulsars
	int *sky_chandra; //Info about if the position of the pulsar was observed by chandra
	int *sky_XMM; //Info about if the position of the pulsar was observed by XMM-Newton
	double *Smin_pmps;
	double *Smin_fast;
	double *w_r_fast;
	double *w_r_pmps;
	double *w_int; 
	double *temp; //Stores sky temperature at (gl,gb) position
	double *Smin_fermi; //Minimum flux detectable by Fermi/LAT
	double **fomega; //f_omega values taken from table
	double *n_omega_x;
	double *n_omega_y;
	double *n_omega_z;
	double *n_mu_x;
	double *n_mu_y;
	double *n_mu_z;
	double *n_nx;
	double *n_ny;
	double *n_nz;
	double *n_sx;
        double *n_sy;
        double *n_sz;
	double *ex;
	double *ey;
	double *ez;
	double *nx;
	double *ny;
	double *nz;
	double *jn; //Angle between rotation axis and <<north>> hotspot
	double *js; //Angle between rotation axis and <<south>> hotspot 
	double *cos_in; //Scalar product of n_in with n_line_of_sight (n)
	double *cos_is; //Scalar product of n_is with n_line_of_sight (n)
	double *theta_n; //saves the localisation of the <<north>> hotspot to check if it is indeed in the north
	double *theta_s; //saves the localisation of the <<south>> hotspot to check if it is indeed in the south
	double *phi_n; //loc <<north>> hotspot
	double *phi_s; //loc <<south>> hotspot
	double *mu_hot_ang1; //angle between <<north>> hotspot and magnetic axis
	double *mu_hot_ang2; //angle between <<south>> hotspot and magnetic axis
	int *see_n; //1 or 0 depending on if we see the <<north>> hotspot
	int *see_s; //1 or 0 depending on if we see the <<south>> hotspot
	double *r_hn; //Radius of the hot spot <<north>>
	double *r_hs; //Radius of the hot spot <<south>>
	double *Temp_n; //Temperature of the hot spot <<north>>
	double *Temp_s; //Temperature of the hot spot <<south>>
	double *PF; //X-ray pulsed fraction
	double *PA;
	double *delta; //Stores the gamma-ray peak separation
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
	double Bfield_const;
	double tau_nu;
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
	double *Fx; //Thermal X-ray flux
	double *Fxmax; //Max of the thermal X-ray flux
	double *Fxmin; //Min of the termal X-ray flux
	double *Fg; //gamma
        double *cos_a0; //a0= initial inclination angle
 	double tau_MHD_al;
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
#define R_NS 10000   /*Moment of Inertia in kg.m2 */ 
#define M_PI 3.14159265358979323846 /* pi */
#define SI_mu0 1.25663706212e-6 /* vacuum permeability in H/m */
#define G_grav 6.67430e-20 //gravitational constant km^-3 kg^-1 s^-2
#define MSUN 1.98847e30 //Solar mass kg
#define SI_eps0 8.85418782e-12 //Vacuum permittivity in F/m
#endif

