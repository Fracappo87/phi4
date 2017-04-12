/* Wilson_Flow: library for flow time integral implementation*/



#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include"Wilson_Flow.h"
#include"geometry.h"
#include"lattice_config.h"





/* WILSON_FLOW_INIT: function for Wilson Flow vector initialization        */
/* -initialize the vectors used for Wilson flow Evolution                  */
/* -initialize Runge-Kutta parameters(Kutta Simpsons' choice in this case) */
/* -initialize the flow equation type                                      */
/* -report all the informations to the logfile                             */



void Wilson_Flow_init(){
    
  Delta_1=(double*)malloc(Vol*sizeof(double));
  Delta_2=(double*)malloc(Vol*sizeof(double));
  Delta_3=(double*)malloc(Vol*sizeof(double));
  Delta_4=(double*)malloc(Vol*sizeof(double));
  Xi_p=(double*)malloc(Vol*sizeof(double));
  
  RK.a=1.0/6.0;
  RK.b=2*RK.a;
  RK.c=RK.b;
  RK.d=RK.a;
  RK.m=0.5;
  RK.lambdar=0.5;
  RK.mu=1.0;
  RK.rho=0.5;
  RK.sigma=1.0;
  RK.tau=0;

  RK.flow_mass=Bulk_mass;
  RK.flow_lambda=Bulk_lambda;
  RK.flow_g=Bulk_g;
  RK.flow_kappa=Bulk_kappa;

  mylog("[WILSON FLOW]: SETTING WILSON FLOW PARAMETERS \n");
  mylog("[WILSON FLOW]: Integration step, dt = %lf\n",dt);
  mylog("[WILSON FLOW]: Number of integration steps, N_steps = %d\n",N_steps);
  mylog("[WILSON FLOW]: Measure rate during flow evolution, N_flow = %d\n",N_flow);
  mylog("[WILSON FLOW]: Wilson Flow implemented with  RK%d integrator \n",int_ord);
  mylog("\n");
  mylog("[WILSON FLOW]: SETTING BULK ACTION PARAMETERS\n");
  mylog("[WILSON FLOW]: Bulk mass parameter, Bulk_mass = %lf\n",RK.flow_mass);
  mylog("[WILSON FLOW]: Bulk lambda parameter, Bulk_lambda = %lf\n",RK.flow_lambda);
  mylog("[WILSON FLOW]: Bulk g parameter, Bulk_g  = %lf\n",RK.flow_g);
  mylog("[WILSON FLOW]: Bulk kappa parameter, Bulk_kappa  = %lf\n",RK.flow_kappa);
  mylog("\n");
  mylog("[WILSON FLOW]: SETTING RUNGE KUTTA PARAMETERS\n");
  mylog("[WILSON FLOW]: Parameter a = %lf\n",RK.a);
  mylog("[WILSON FLOW]: Parameter b = %lf\n",RK.b);
  mylog("[WILSON FLOW]: Parameter c = %lf\n",RK.c);
  mylog("[WILSON FLOW]: Parameter d = %lf\n",RK.d);
  mylog("[WILSON FLOW]: Parameter m = %lf\n",RK.m);
  mylog("[WILSON FLOW]: Parameter lambda = %lf\n",RK.lambdar);
  mylog("[WILSON FLOW]: Parameter mu = %lf\n",RK.mu);
  mylog("[WILSON FLOW]: Parameter rho = %lf\n",RK.rho);
  mylog("[WILSON FLOW]: Parameter sigma = %lf\n",RK.sigma);
  mylog("[WILSON FLOW]: Parameter tau = %lf\n",RK.tau);

  mylog("\n");
  
}





/* FLOW_EQ:f unction used to evaluate the discretized flow equation at every time slice */



double Flow_eq(double *vec,int i){

  int mu;
  double flow_eq;

  flow_eq=-RK.flow_mass*RK.flow_mass*vec[i]-(RK.flow_g/2.0)*vec[i]*vec[i]-(RK.flow_lambda/6.0)*vec[i]*vec[i]*vec[i]-(RK.flow_kappa/120.0)*vec[i]*vec[i]*vec[i]*vec[i]*vec[i];

  for(mu=0;mu<d;mu++){
    
    flow_eq+=vec[d_up[i*d+mu]]+vec[d_down[i*d+mu]]-2*vec[i];  

  }

  return flow_eq;

}





/* WILSON_FLOW_RK1: 1st order Runge-Kutta algorithm for the integration of the flow equation                                                    */
/* -define 1st order function increment, time evolved field Xi                                                                                  */
/* Single step at a given slice time t:                                                                                                         */
/*  L-for every lattice site i                                                                                                                  */
/*  -define Delta_1[i]=Flow_eq(Xi,i)*dt                                                                                                         */
/*  -define Xi[i]=Xi[i]+Delta_1[i]                                                                                                              */



void Wilson_Flow_RK1(double delta_t, double *Xi){
  
  int i;
  
  for(i=0;i<Vol;i++){
    
    Delta_1[i]=Flow_eq(Xi,i)*delta_t;
    Xi[i]=Xi[i]+Delta_1[i];
    
  }
}





/* WILSON_FLOW_RK2: 2nd order Runge-Kutta algorithm for the integration of the flow equation                                                               */
/* -define 1st, 2nd order function increments, time evolved field Xi, an auxialiry field Xi_p used for increment evaluation                                */
/* Single step at a given slice time t:                                                                                                                    */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_1[i]=Flow_eq(Xi,i)*dt                                                                                                                    */
/*  -define  Xi_p[i]=Xi[i]+Delta_1[i]( we need this vector for the evaluation of Delta_2)                                                                  */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_2[i]=Flow_eq(Xi_p,i)*dt                                                                                                                  */
/*  -define Xi[i]=Xi[i]+0.5*(Delta_1[i]+Delta_2[i]                                                                                                         */



void Wilson_Flow_RK2(double delta_t, double *Xi){
  
  int i;
  
  for(i=0;i<Vol;i++){
    
    Delta_1[i]=Flow_eq(Xi,i)*delta_t;
    Xi_p[i]=Xi[i]+Delta_1[i];

  }
  
   for(i=0;i<Vol;i++){
    
    Delta_2[i]=Flow_eq(Xi_p,i)*delta_t;
    Xi[i]=Xi[i]+0.5*(Delta_1[i]+Delta_2[i]);
    
    }
}





/* WILSON_FLOW_RK3: 3rd order Runge-Kutta algorithm for the integration of the flow equation                                                               */
/* -define 1st, 2nd,3rd order function increments, time evolved field Xi, an auxialiry field Xi_p used for increment evaluation                            */
/* Single step at a given slice time t:                                                                                                                    */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_1[i]=Flow_eq(Xi,i)*dt                                                                                                                    */
/*  -define  Xi_p[i]=Xi[i]+0.5*Delta_1[i]( we need this vector for the evaluation of Delta_2)                                                              */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_2[i]=Flow_eq(Xi_p,i)*dt                                                                                                                  */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Xi_p[i]=Xi[i]-1.0*Delta_1[i]+2.0*Delta_2[i]( we need this vector for the evaluation of Delta_3)                                                */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_3[i]=Flow_eq(Xi_p,i)*dt                                                                                                                  */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Xi[i]=Xi[i]+(1.0/6.0)*(Delta_1[i]+4.0*Delta_2[i]+Delta_3[i]                                                                                    */



void Wilson_Flow_RK3(double delta_t, double *Xi){
  
  int i;
  
  for(i=0;i<Vol;i++){
    
    Delta_1[i]=Flow_eq(Xi,i)*delta_t;
    Xi_p[i]=Xi[i]+0.5*Delta_1[i];
    
  }
  
  for(i=0;i<Vol;i++) Delta_2[i]=Flow_eq(Xi_p,i)*delta_t;
  for(i=0;i<Vol;i++) Xi_p[i]=Xi[i]-1.0*Delta_1[i]+2.0*(Delta_2[i]);

  for(i=0;i<Vol;i++) Delta_3[i]=Flow_eq(Xi_p,i)*delta_t;
  for(i=0;i<Vol;i++) Xi[i]+=(1.0/6.0)*(Delta_1[i]+4.0*Delta_2[i]+Delta_3[i]);
  
}





/* WILSON_FLOW_RK4: 4th order Runge-Kutta algorithm for the integration of the flow equation                                                               */
/* -define 1st, 2nd...4th order function increments, time evolved field Xi, an auxialiry field Xi_p used for increment evaluation                          */
/* Single step at a given slice time t:                                                                                                                    */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_1[i]=Flow_eq(Xi,i)*dt                                                                                                                    */
/*  -define  Xi_p[i]=Xi[i]+(RK.m)*Delta_1[i]( we need this vector for the evaluation of Delta_2)                                                           */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_2[i]=Flow_eq(Xi_p,i)*dt                                                                                                                  */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Xi_p[i]=Xi[i]+(RK.lambda-RK.rho)*Delta_1[i]+(RK.rho)*Delta_2[i]( we need this vector for the evaluation of Delta_3)                            */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_3[i]=Flow_eq(Xi_p,i)*dt                                                                                                                  */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Xi_p[i]=Xi[i]+(RK.mu-RK.sigma-RK.tau)*Delta_1[i]+(RK.tau)*Delta_2[i]+(RK.sigma)*Delta_3[i]( we need this vector for the evaluation of Delta_4) */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Delta_4[i]=Flow_eq(Xi_p,i)*dt                                                                                                                  */
/*  L-for every lattice site i                                                                                                                             */
/*  -define Xi[i]=Xi[i]+(RK.a)*Delta_1[i]+(RK.b)*Delta_2[i]+(RK.c)*Delta_3[i]+(RK.d)*Delta_4[i]                                                            */
/* -different choiches of the coefficents give a different Runge Kutta integrator: they can be, CAREFULLY, modified in Wilson_Flow_init                    */



void Wilson_Flow_RK4(double delta_t,double *Xi){
  
  int i;

    for(i=0;i<Vol;i++){
      
      Delta_1[i]=Flow_eq(Xi,i)*delta_t;
      Xi_p[i]=Xi[i]+(RK.m)*Delta_1[i];
      
    }
    
    for(i=0;i<Vol;i++) Delta_2[i]=Flow_eq(Xi_p,i)*delta_t;
    for(i=0;i<Vol;i++) Xi_p[i]=Xi[i]+(RK.lambdar-RK.rho)*Delta_1[i]+(RK.rho)*Delta_2[i];
    
    for(i=0;i<Vol;i++) Delta_3[i]=Flow_eq(Xi_p,i)*delta_t;
    for(i=0;i<Vol;i++) Xi_p[i]=Xi[i]+(RK.mu-RK.sigma-RK.tau)*Delta_1[i]+(RK.tau)*Delta_2[i]+(RK.sigma)*Delta_3[i];
    
    for(i=0;i<Vol;i++) Delta_4[i]=Flow_eq(Xi_p,i)*delta_t;
    for(i=0;i<Vol;i++) Xi[i]+=(RK.a)*Delta_1[i]+(RK.b)*Delta_2[i]+(RK.c)*Delta_3[i]+(RK.d)*Delta_4[i];
    
}





/* FREE_WF: function used for freeing Wilson Flow vectors memory */



void free_WF(double *Xi){

  free(Delta_1);
  Delta_1=NULL;
  free(Delta_2);
  Delta_2=NULL;
  free(Delta_3);
  Delta_3=NULL;
  free(Delta_4);
  Delta_4=NULL;
  free(Xi_p);
  Xi_p=NULL;
  free(Xi);
  Xi=NULL;

  mylog("[WILSON FLOW]: Auxiliary fields' memory freed\n");
  mylog("[WILSON FLOW]: Flow evolved field's memory freed\n");
  mylog("\n");

}
