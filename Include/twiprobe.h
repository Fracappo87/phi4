
#ifndef __TWIPROB
#define __TWIPROB

#include <math.h>
#include "macro.h"
#include <stdio.h>
#include "observables.h"
#include "geometry.h"
#include "lattice_config.h"


#define NUM_PROB 11

typedef double (*prb_ptr_t)(int,int,double *);
double probe_1(int x, int alpha,double *Phi);
double probe_2(int x, int alpha,double *Phi);
double probe_3(int x, int alpha,double *Phi);
double probe_4(int x, int alpha,double *Phi);
double probe_5(int x, int alpha,double *Phi);
double probe_6(int x, int alpha,double *Phi);
double probe_7(int x, int alpha,double *Phi);
double probe_8(int x, int alpha,double *Phi);
double probe_9(int x, int alpha,double *Phi);
double probe_10(int x, int alpha,double *Phi);
double probe_11(int x, int alpha,double *Phi);

#endif