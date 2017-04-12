/* OBSERVABLES: library for the definitions of measured observables */



#ifndef __OBSERVABLES
#define __OBSERVABLES

#include <math.h>
#include "lattice_config.h"
#include "complex.h"
#include "Jacobian.h"
#define POWER 2



/* Obs: macro we use to define the probe observable of the integrated TWI*/


#ifndef POWER
#define POWER 1
#endif
#if POWER==1
#define Obs(site)  phi[(site)] 
#define DObs(site) 1.
#endif
#if POWER==2
#define Obs(site)  phi[(site)]*phi[(site)] 
#define DObs(site) 2.*phi[(site)]
#endif

/* MEAN: define the mean */



double mean(double* vec, int d);





/* CORR: function for the evaluation of correlators */



double corr(int dis,double *phi);
double rhschainnew(int dis, double *phi);




/* SPATIAL_MEAN: d-1 dimensional mean, used for TWI */



double spatial_mean(int coord_0, double *phi);





/* SPATIAL_MEAN_MUL: d-1 dimensional mean of operators product, used for TWI */



double spatial_mean_mul(int coord_0,double *phi,double *psi);





/* SPATIAL_MEAN_MUL2: d-1 dimensional mean of operators product, used for TWI */



double spatial_mean_mul2(int coord_0,double *phi);





/* SPATIAL_MEAN: d-1 dimensional mean of operators product, used for TWI */



double spatial_mean_mul3(int coord_0, double *phi);





/* PHI2: define the mean of a squared lattice quantity */



double phi2(double *phi);





/* PHI4: define the mean of the fourth power of a lattice quantity */



double phi4(double *phi);





/* PHI6: define the mean of a sixth power of a lattice quantity */



double phi6(double *phi);





/* PHI_TILDA: define the fourier transform of the field */



void phi_tilda(double* p, double* phi,complex* phitilda);





/* ACTION: define the action of the system */



double action(double *phi);





/* EMT: function for the evalutation of the EMT components */



double EMT(int mu, int nu, int n, double *Field);





/* T1: mixing operator with Tmunu */



double T1(int i, int nu);





/* T2: mixing operator with Tmunu */



double T2(int i, int nu);





/* T3: mixing operator with Tmunu */



double T3(int i, int nu);





/* T4: mixing operator with Tmunu */



double T4(int i, int nu);





/* T5: mixing operator with Tmunu */



double T5(int i, int nu,int rho);





/* ZRHS: r.h.s of the integrated TWI with naive lattice derivative*/



double Zrhs(int i0,int j0,double *phi,char name);





/* ZRHS: r.h.s of the integrated TWI with chain rule lattice derivative */



double Zrhs_chain(int i0,int j0,double *phi, double *dphi);





/* VIOLATION: function that measure the numerical difference from naive and chain rule lattice derivative formulations*/



double violation(int i0, int j0, double *phi,char name);





/* ZLHS: l.h.s. of the integrated TWI */



double Zlhs(int i0, int j0, double *phi, double *psi);





/* LHPERT: alterinative version of the l.h.s. of the integrated TWI */



double LHpert(int coord_0,int coord_1, double *phi, double *psi);





/* LHPERT_exact: alterinative version of the l.h.s. of the integrated TWI, using the exact form of the Jacobian */



double LHpert_exact(int coord_0,int coord_1,int size,int t, double *dphi, double *phi, double *jac_vec);
double LHpert_exact_vol(int dis,int size,int t, double *dphi, double *phi, double *jac_vec);

#endif
