/* Wilson_Flow: library for flow time integral implementation*/



#ifndef __WILSON_FLOW
#define __WILSON_FLOW

#include"glob_const.h"




/* DELTAs, XIs: vector used for Wilson flow implementation */



double *Delta_1;
double *Delta_2;
double *Delta_3;
double *Delta_4;
double *Xi_p;






/* WILSON_FLOW_PARAMETERS: structure used to include the value of Runge Kutta  parameters and define the flow equation type*/



typedef struct Wilson_Flow_Parameters{

  double a;
  double b;
  double c;
  double d;
  double m;
  double lambdar;
  double mu;  
  double rho;
  double sigma;
  double tau;
 
  double flow_mass;
  double flow_lambda;
  double flow_g;
  double flow_kappa;

}WFP;

WFP RK;





/* WILSON_FLOW_INIT: function for Wilson Flow vector and Runge Kutta parameters initialization */



void Wilson_Flow_init();





/* FLOW_EQ:f unction used to evaluate the discretize flow equation at every time slice */



double Flow_eq(double *vec, int i);





/* WILSON_FLOW_RK1: 1st order Runge-Kutta algorithm for the integration of the flow equation */



void Wilson_Flow_RK1(double delta_t, double *Xi);





/*WILSON_FLOW_2ND:  2nd order Runge-Kutta algorithm for the integration of the flow equation */



void Wilson_Flow_RK2(double delta_t, double *Xi);



  

/*WILSON_FLOW_RK3:  3rd order Runge-Kutta algorithm for the integration of the flow equation */



void Wilson_Flow_RK3(double delta_t, double *Xi);





/* WILSON_FLOW_RK4: 4th order Runge-Kutta algorithm for the integration of the flow equation */



void Wilson_Flow_RK4(double delta_t, double *Xi);





/* FREE_WF: function used for freeing Wilson Flow vectors memory */



void free_WF(double *Xi);


#endif
