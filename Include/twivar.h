
#ifndef __TWIVAR
#define __TWIVAR

#include <math.h>
#include "macro.h"
#include "macro_jacobian.h"
#include <stdio.h>
#include "observables.h"
#include "geometry.h"
#include "lattice_config.h"


typedef double (*dprb_ptr_t)(int,int,int,int,int,double *, double *);

double TWIJac(int x,int y,int tflow);
double delta_pr_1(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_2(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_3(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_4(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_5(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_6(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_7(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_8(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_9(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_10(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);
double delta_pr_11(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary);

#endif