/* Jacobian: library for constructing the Jacobian matrix related to the flow equation */



#ifndef __JACOBIAN
#define __JACOBIAN

#include"glob_const.h"
#include"Wilson_Flow.h"




/* JACOBIAN_FLOW_EQ:f unction used to evaluate the discretize flow equation of the Jacobian at every time slice */



double Jacobian_Flow_eq(double *vec_1, double *vec_2, int i);





/* ZDELTA_SET_CONDITION: function used to employ the action of the Jacobian operator, reduced to a subdomain D of the global lattive volume V */



void Zdelta_set_condition(int coord_0,int d,double *J);





/* JACOBIAN_FLOW_RK1: 1st order Runge-Kutta algorithm for the integration of the Jacobian flow equation */



void Jacobian_Flow_RK1(double delta_t, double *Xi, double *J);





/*JACOBIAN_FLOW_2ND: 2nd order Runge-Kutta algorithm for the integration of the Jacobian flow equation */



void Jacobian_Flow_RK2(double delta_t, double *Xi, double *J);



  

/*JACOBIAN_FLOW_RK3: 3rd order Runge-Kutta algorithm for the integration of the Jacobian flow equation  */



void Jacobian_Flow_RK3(double delta_t, double *Xi, double *J);





/* JACOBIAN_FLOW_RK4: 4th order Runge-Kutta algorithm for the integration of the Jacobian flow equation */



void Jacobian_Flow_RK4(double delta_t, double *Xi, double *J);





/* INIT_JAC: numerical tabulation of the exact lattice Jacobian (other possible way to implement the TWI in the pure smearing case)*/



double *jac_v;
double **Jacobian;
int DIST;
void jac(int x1,int x2,int x3);
void jac_rk(int x1,int x2,int x3,double dt);
void jac_onepoint(int x1,int x2,int x3);


double *dt_jac_v;
void dt_jac(void);
void dt_jac_onepoint(int x1,int x2,int x3);



/* FREE_J: function used for freeing Jacobian vector memory */



void free_J(double *J);





/* FREE_JAC: function for freeing the Jacobian table  */



void free_jac(void);

void free_dt_jac(void);





#endif
