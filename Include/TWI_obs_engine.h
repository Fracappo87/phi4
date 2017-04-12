/* TWI_OBS_ENGINE: routine for the definitions of all the operators that enters in the definition of our local translation Ward Identities*/

#ifndef __TWI_OBS
#define __TWI_OBS

#include <math.h>
#include "macro.h"
#include <stdio.h>
#include "observables.h"
#include "geometry.h"
#include "lattice_config.h"
#include "twiinclude.h"
#define NUM_OP 4




typedef double (*obs_ptr_t)(int,int,int,double *);


double *avg_vector_init(void);
void vec_clearer(double *vec);

double phi_1(int i, int mu, int nu, double *phi);
double phi_2(int i, int mu, int nu, double *phi);
double phi_3(int i, int mu, int nu, double *phi);
double phi_4(int i, int mu, int nu, double *phi);
double phi_5(int i, int mu, int nu, double *phi);
double phi_6(int i, int mu, int nu, double *phi);
double phi_7(int i, int mu, int nu, double *phi);
double phi_8(int i, int mu, int nu, double *phi);

double LHS_delta_V(int nprob,int x,int y, int tflow, double *phiflow,double *phi);

double RHS_loc(int x,int y,int nprob,int ntens,double *phiflow,double *phi);
double RHS(int x,int y,int nprob,int ntens,double *phiflow,double *phi);


#endif
