/* JACOBTEST: library used for testing Runge Kutta integrators for the Jacobian matrix */



#ifndef __JACOBTEST
#define __JACOBTEST

#include"glob_const.h"
#include"testRK.h"





/* DRHO_PHIANALYTICCOS: function used to define analytic solution to be compared with the numerical one */



double drho_phianalyticcos(double t, int i,int rho, double *p);





/* TESTJACOB_INIT: function used to set the initial condition of the numerical solution. */



void testJACOB_init(double *p, double *X, double *J);





/* TESTJACOB: function used to test the efficency of a given integration algorithm */



void testJACOB();





#endif
