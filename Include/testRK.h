/* RKTEST: library used for testing Runge Kutta integrators */



#ifndef __RKTEST
#define __RKTEST

#include"glob_const.h"





/* LATTICE_MOM: trivial function we need in order to define the lattice momentum (see any lattice field theory reference for more infos )*/



double lattice_mom(double *p,int mu);





/* SQ_LATTICE_MOM: trivial function we need in order to define the square lattice momentum (see any lattice field theory reference for more infos )*/



double sq_lattice_mom(double *p);





/* PHIANALYTICCOS: function used to define the analytic solution to be compared with the numerical one */



double phianalyticcos(double t, int i, double *p);





/* TESTRK_INIT: function used to set the initial condition of the numerical solution. */



void testRK_init(double *p,double *X);





/* TESTRK: function used to test the efficency of a given integration algorithm */



void testRK();


#endif
