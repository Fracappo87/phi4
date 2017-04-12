/* DWI_OBS_ENGINE: routine for the definitions of all the operators that enters in the definition of our local translation Ward Identities*/

#ifndef __DWI_OBS
#define __DWI_OBS

#include <math.h>
#include <stdio.h>
#include "macro.h"
#include "observables.h"
#include "geometry.h"
#include "lattice_config.h"

#define dwiNUM_PROB 8
#define dwiNUM_OP 5

typedef double (*dwiobs_ptr_t)(int,double *);
typedef double (*dwiprb_ptr_t)(int,double *);
typedef double (*dwiprb_der_ptr_t)(int,int,double*);

double *dwiavg_vector_init(void);
void dwivec_clearer(double *vec);

double dwiphi_2(int i,  double *phi);
double dwiphi_4(int i, double *phi);
double dwiphi_5(int i,  double *phi);
double dwiphi_6(int i,  double *phi);
double dwiphi_0(int i,  double *phi);

double dwipr_1(int i,double *phi);
double dwipr_2(int i,double *phi);
double dwipr_3(int i,double *phi);
double dwipr_4(int i,double *phi);
double dwipr_5(int i,double *phi);
double dwipr_6(int i,double *phi);
double dwipr_7(int i,double *phi);
double dwipr_8(int i,double *phi);

double dphi_pr_1(int i,double *phi);
double dphi_pr_2(int i,double *phi);
double dphi_pr_3(int i,double *phi);
double dphi_pr_4(int i,double *phi);
double dphi_pr_5(int i,double *phi);
double dphi_pr_6(int i,double *phi);
double dphi_pr_7(int i,double *phi);
double dphi_pr_8(int i,double *phi);

double ddmuphi_pr_1(int i,int rho,double *phi);
double ddmuphi_pr_2(int i,int rho,double *phi);
double ddmuphi_pr_3(int i,int rho,double *phi);
double ddmuphi_pr_4(int i,int rho,double *phi);
double ddmuphi_pr_5(int i,int rho,double *phi);
double ddmuphi_pr_6(int i,int rho,double *phi);
double ddmuphi_pr_7(int i,int rho,double *phi);
double ddmuphi_pr_8(int i,int rho,double *phi);

double dt_phi(int i, int tflow, double *phi);
double dtdmu_phi(int i,int mu,int tflow,double *phi);

double dwiLHS(int npro,int y, int tflow, double *phi_flow,double *phi);
double dwiRHS(int nprob,int ntens,int y,double *phi_flow,double *phi);

#endif
