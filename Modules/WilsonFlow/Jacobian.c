/* Jacobian: library for constructing the Jacobian matrix related to the flow equation */



#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include "Wilson_Flow.h"
#include"Jacobian.h"
#include"geometry.h"
#include"lattice_config.h"





/* JACOBIAN_FLOW_EQ:f unction used to evaluate the discretized flow equation of the Jacobian at every time slice */



double Jacobian_Flow_eq(double *vec_1, double *vec_2, int i){

  int mu;
  double flow_eq;

  flow_eq=(-RK.flow_mass*RK.flow_mass-RK.flow_g*vec_1[i]-(RK.flow_lambda/2.0)*vec_1[i]*vec_1[i]-(RK.flow_kappa/24.0)*vec_1[i]*vec_1[i]*vec_1[i]*vec_1[i])*vec_2[i];

  for(mu=0;mu<d;mu++){
    
    flow_eq+=vec_2[d_up[i*d+mu]]+vec_2[d_down[i*d+mu]]-2*vec_2[i];  

  }
  
  return flow_eq;

}





/* Zdelta_set_condition: this function forces the action of the Jacobian, along its all time evolution, to be restricted to a subdomain of the entire volume */
/* The easiest way to implement such a constriction is customizing the initial condition on the integral action of the Jacobian                              */
/* In this case, our domani is given by [x0-l,x0+l]x[0,L]^{d-1}, where L is the number of lattice point in a given direction, l <=L                          */
/* More precise informations can be found in the article of Del Debbio, Patella, Rago [http://arxiv.org/abs/1306.1173]                                       */



void Zdelta_set_condition(int coord_0,int d,double *J){

  int i;
  int x[d];
  
  for(i=0;i<Vol;i++){
    
    coordinates(i,x);
    
    if(abs(coord_0-x[0]) > d && abs(-coord_0+x[0]+Nt) >d) J[i]=0.;
  }
  
}





/* JACOBIAN_FLOW_RK1: 1st order Runge-Kutta algorithm for the integration of the Jacobian flow equation                                         */
/* -define 1st order function increment, time evolved field Xi                                                                                  */
/* Single step at a given slice time t:                                                                                                         */
/*  L-for every lattice site i                                                                                                                  */
/*  -define Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*dt                                                                                              */
/*  -define J[i]=J[i]+Delta_1[i]                                                                                                                */



void Jacobian_Flow_RK1(double delta_t, double *Xi, double *J){
  
  int i;
  
  for(i=0;i<Vol;i++){
    
    Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*delta_t;
    J[i]=J[i]+Delta_1[i];
    
  }
}





/* JACOBIAN_FLOW_RK2: 2nd order Runge-Kutta algorithm for the integration of the Jacobian flow equation                                                   */
/* -define 1st, 2nd order function increments, time evolved field J, an auxialiry field Xi_p used for increment evaluation                                */
/* Single step at a given slice time t:                                                                                                                   */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*dt                                                                                                        */
/*  -define  Xi_p[i]=J[i]+Delta_1[i]( we need this vector for the evaluation of Delta_2)                                                                  */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_2[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*dt                                                                                                     */
/*  -define J[i]=J[i]+0.5*(Delta_1[i]+Delta_2[i]                                                                                                          */



void Jacobian_Flow_RK2(double delta_t, double *Xi, double *J){
  
  int i;
  
  for(i=0;i<Vol;i++){
    
    Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*delta_t;
    Xi_p[i]=J[i]+Delta_1[i];

  }
  
   for(i=0;i<Vol;i++){
    
    Delta_2[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*delta_t;
    J[i]=J[i]+0.5*(Delta_1[i]+Delta_2[i]);
    
    }
}





/* JACOBIAN_FLOW_RK3: 3rd order Runge-Kutta algorithm for the integration of the Jacobian flow equation                                                   */
/* -define 1st, 2nd,3rd order function increments, time evolved field J, an auxialiry field Xi_p used for increment evaluation                            */
/* Single step at a given slice time t:                                                                                                                   */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*dt                                                                                                        */
/*  -define  Xi_p[i]=J[i]+0.5*Delta_1[i]( we need this vector for the evaluation of Delta_2)                                                              */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_2[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*dt                                                                                                     */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Xi_p[i]=J[i]-1.0*Delta_1[i]+2.0*Delta_2[i]( we need this vector for the evaluation of Delta_3)                                                */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_3[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*dt                                                                                                     */
/*  L-for every lattice site i                                                                                                                            */
/*  -define J[i]=J[i]+(1.0/6.0)*(Delta_1[i]+4.0*Delta_2[i]+Delta_3[i]                                                                                     */



void Jacobian_Flow_RK3(double delta_t, double *Xi, double *J){
  
  int i;
  
  for(i=0;i<Vol;i++){
    
    Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*delta_t;
    Xi_p[i]=J[i]+0.5*Delta_1[i];
    
  }
  
  for(i=0;i<Vol;i++) Delta_2[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*delta_t;
  for(i=0;i<Vol;i++) Xi_p[i]=J[i]-1.0*Delta_1[i]+2.0*(Delta_2[i]);

  for(i=0;i<Vol;i++) Delta_3[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*delta_t;
  for(i=0;i<Vol;i++) J[i]+=(1.0/6.0)*(Delta_1[i]+4.0*Delta_2[i]+Delta_3[i]);
  
}





/* JACOBIAN_FLOW_RK4: 4th order Runge-Kutta algorithm for the integration of the Jacobian flow equation                                                   */
/* -define 1st, 2nd...4th order function increments, time evolved field J, an auxialiry field Xi_p used for increment evaluation                          */
/* Single step at a given slice time t:                                                                                                                   */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*dt                                                                                                        */
/*  -define  Xi_p[i]=J[i]+(RK.m)*Delta_1[i]( we need this vector for the evaluation of Delta_2)                                                           */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_2[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*dt                                                                                                     */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Xi_p[i]=J[i]+(RK.lambda-RK.rho)*Delta_1[i]+(RK.rho)*Delta_2[i]( we need this vector for the evaluation of Delta_3)                            */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_3[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*dt                                                                                                     */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Xi_p[i]=J[i]+(RK.mu-RK.sigma-RK.tau)*Delta_1[i]+(RK.tau)*Delta_2[i]+(RK.sigma)*Delta_3[i]( we need this vector for the evaluation of Delta_4) */
/*  L-for every lattice site i                                                                                                                            */
/*  -define Delta_4[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*dt                                                                                                     */
/*  L-for every lattice site i                                                                                                                            */
/*  -define J[i]=J[i]+(RK.a)*Delta_1[i]+(RK.b)*Delta_2[i]+(RK.c)*Delta_3[i]+(RK.d)*Delta_4[i]                                                             */
/* -different choiches of the coefficents give a different Runge Kutta integrator: they can be, CAREFULLY, modified in Wilson_Flow_init                   */



void Jacobian_Flow_RK4(double delta_t,double *Xi, double *J){
  
  int i;

    for(i=0;i<Vol;i++){
      
      Delta_1[i]=Jacobian_Flow_eq(Xi,J,i)*delta_t;
      Xi_p[i]=J[i]+(RK.m)*Delta_1[i];
      
    }
    
    for(i=0;i<Vol;i++) Delta_2[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*delta_t;
    for(i=0;i<Vol;i++) Xi_p[i]=J[i]+(RK.lambdar-RK.rho)*Delta_1[i]+(RK.rho)*Delta_2[i];
    
    for(i=0;i<Vol;i++) Delta_3[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*delta_t;
    for(i=0;i<Vol;i++) Xi_p[i]=J[i]+(RK.mu-RK.sigma-RK.tau)*Delta_1[i]+(RK.tau)*Delta_2[i]+(RK.sigma)*Delta_3[i];
    
    for(i=0;i<Vol;i++) Delta_4[i]=Jacobian_Flow_eq(Xi,Xi_p,i)*delta_t;
    for(i=0;i<Vol;i++) J[i]+=(RK.a)*Delta_1[i]+(RK.b)*Delta_2[i]+(RK.c)*Delta_3[i]+(RK.d)*Delta_4[i];
    
}





/* JAC: numerical tabulation of the exact lattice Jacobian (other possible way to implement the TWI in the pure smearing case)*/



void jac(int x1,int x2,int x3){

  int i;
  double *mom,*latmom; 
  DIST=get_index(x1,x2,x3);
  mom=(double*) malloc(d*Vol*sizeof(double));
  latmom=(double*) malloc(d*Vol*sizeof(double));
  jac_v=(double*) malloc(Vol*(N_steps+1)*sizeof(double));  
  Jacobian=(double **) malloc((N_steps+1)*sizeof(double *));
  mylog("[JACOBIAN]: Starting Jacobian tabulation \n");
  mylog("\n");

  int t_fl,index,iext;
  for(t_fl=0;t_fl<=N_steps;t_fl++){
    Jacobian[t_fl]=jac_v + t_fl * Vol;

  }
  
  
  for(i=0;i<Vol;i++){
    int x[d];
    coordinates(i,x);
    mom[i]=2.*mPI*x[0]/(double)Nt;
    mom[i+Vol]=2.*mPI*x[1]/(double)Nx;
    mom[i+2*Vol]=2.*mPI*x[2]/(double)Ny;
    
    latmom[i]=2.*sin(mom[i]/2.);
    latmom[i+Vol]=2.*sin(mom[i+Vol]/2.);
    latmom[i+2*Vol]=2.*sin(mom[i+2*Vol]/2.);
    
  }
  
  

  double sum;
  double scalarprod;
  double latmomsq;
  
  for(iext=0;iext<Vol;iext++){
    int x[d];
    coordinates(iext,x);
    for(t_fl=0;t_fl<=N_steps;t_fl++){
      
      sum=0.;
      
      for(index=0;index<Vol;index++){
	scalarprod=mom[index]*x[0]+mom[index+Vol]*x[1]+mom[index+2*Vol]*x[2];
	latmomsq=latmom[index]*latmom[index]+latmom[index+Vol]*latmom[index+Vol]+latmom[index+2*Vol]*latmom[index+2*Vol];
	sum+=cos(scalarprod)*exp(-t_fl*dt*latmomsq);
      }
      
      sum/=(double)Vol;
      jac_v[iext+Vol*t_fl]=sum;
      
    }
  }
    
  free(latmom);
  free(mom);
  
  mylog("[JACOBIAN]: Tabulation ended \n");
  mylog("\n");
  
    return;
}



void jac_rk(int x1,int x2,int x3,double dt){

  DIST=get_index(x1,x2,x3);
  jac_v=(double*) malloc(Vol*(N_flowtimes+1)*sizeof(double));  
  double *xi=(double *) malloc(Vol*sizeof(double));
  Jacobian=(double **) malloc((N_flowtimes+1)*sizeof(double *));
  mylog("[JACOBIAN]: Starting Jacobian tabulation \n");
  mylog("\n");

  int t_fl,iext;
  for(t_fl=0;t_fl<=N_flowtimes;t_fl++)
    Jacobian[t_fl]=jac_v + t_fl * Vol;
  
  xi[0]=1.;
  jac_v[0]=1.;
  for(iext=1;iext<Vol;iext++){
    xi[iext]=0.;
    jac_v[iext]=0.;
      }

  int count=0;

  for(t_fl=0;t_fl<N_steps;t_fl++){
    Wilson_Flow_RK4(dt,xi);
    if((Flow_vec[count]-t_fl-1)==0){
      for(iext=0;iext<Vol;iext++)
	jac_v[iext+(count+1)*Vol]=xi[iext];
      count++;
    }
  }





  
  
  mylog("[JACOBIAN]: Tabulation ended \n");
  mylog("\n");
  
    return;
}



void jac_onepoint(int x1,int x2,int x3){

  int i;
  double *mom,*latmom; 

  mom=(double*) malloc(d*Vol*sizeof(double));
  latmom=(double*) malloc(d*Vol*sizeof(double));
  jac_v=(double*) malloc((N_steps+1)*sizeof(double));  
  
  mylog("[JACOBIAN]: Starting Jacobian tabulation \n");
  mylog("\n");
  
  
  for(i=0;i<Vol;i++){
    int x[d];
    coordinates(i,x);
    mom[i]=2.*mPI*x[0]/(double)Nt;
    mom[i+Vol]=2.*mPI*x[1]/(double)Nx;
    mom[i+2*Vol]=2.*mPI*x[2]/(double)Ny;
    
    latmom[i]=2.*sin(mom[i]/2.);
    latmom[i+Vol]=2.*sin(mom[i+Vol]/2.);
    latmom[i+2*Vol]=2.*sin(mom[i+2*Vol]/2.);
    
  }
  
  
  int t_fl,index;
  double sum;
  double scalarprod;
  double latmomsq;
  
    for(t_fl=0;t_fl<=N_steps;t_fl++){
      
      sum=0.;
      
      for(index=0;index<Vol;index++){
	scalarprod=mom[index]*x1+mom[index+Vol]*x2+mom[index+2*Vol]*x3;
	latmomsq=latmom[index]*latmom[index]+latmom[index+Vol]*latmom[index+Vol]+latmom[index+2*Vol]*latmom[index+2*Vol];
	sum+=cos(scalarprod)*exp(-t_fl*dt*latmomsq);
      }
      
      sum/=(double)Vol;
      jac_v[t_fl]=sum;
      
    }
    
  free(latmom);
  free(mom);
  
  mylog("[JACOBIAN]: Tabulation ended \n");
  mylog("\n");
  
    return;
}


/* JAC: numerical tabulation of the exact lattice Jacobian (other possible way to implement the TWI in the pure smearing case)*/



void dt_jac(){

  int i;
  double *mom,*latmom; 

  mom=(double*) malloc(d*Vol*sizeof(double));
  latmom=(double*) malloc(d*Vol*sizeof(double));
  dt_jac_v=(double*) malloc(Vol*(N_steps+1)*sizeof(double));  
  



  mylog("[JACOBIAN]: Starting flow derivative of the Jacobian tabulation \n");
  mylog("\n");
  
  
  for(i=0;i<Vol;i++){
    int x[d];
    coordinates(i,x);
    mom[i]=2.*mPI*x[0]/(double)Nt;
    mom[i+Vol]=2.*mPI*x[1]/(double)Nx;
    mom[i+2*Vol]=2.*mPI*x[2]/(double)Ny;
    
    latmom[i]=2.*sin(mom[i]/2.);
    latmom[i+Vol]=2.*sin(mom[i+Vol]/2.);
    latmom[i+2*Vol]=2.*sin(mom[i+2*Vol]/2.);
    
  }
  
  
  int t_fl,index,iext;
  double sum;
  double scalarprod;
  double latmomsq;
  
  for(iext=0;iext<Vol;iext++){
    int x[d];
    coordinates(iext,x);
    for(t_fl=0;t_fl<=N_steps;t_fl++){
      
      sum=0.;
      
      for(index=0;index<Vol;index++){
	scalarprod=mom[index]*x[0]+mom[index+Vol]*x[1]+mom[index+2*Vol]*x[2];
	latmomsq=latmom[index]*latmom[index]+latmom[index+Vol]*latmom[index+Vol]+latmom[index+2*Vol]*latmom[index+2*Vol];
	sum+=-latmomsq*cos(scalarprod)*exp(-t_fl*dt*latmomsq);
      }
      
      sum/=(double)Vol;
      dt_jac_v[iext+Vol*t_fl]=sum;
      
    }
  }
    
  free(latmom);
  free(mom);
  
  mylog("[JACOBIAN]: Tabulation ended \n");
  mylog("\n");
  
    return;
}




void dt_jac_onepoint(int x1,int x2,int x3){

  int i;
  double *mom,*latmom; 

  mom=(double*) malloc(d*Vol*sizeof(double));
  latmom=(double*) malloc(d*Vol*sizeof(double));
  dt_jac_v=(double*) malloc((N_steps+1)*sizeof(double));  
  
  mylog("[JACOBIAN]: Starting flow derivative of the Jacobian tabulation \n");
  mylog("\n");
  
  
  for(i=0;i<Vol;i++){
    int x[d];
    coordinates(i,x);
    mom[i]=2.*mPI*x[0]/(double)Nt;
    mom[i+Vol]=2.*mPI*x[1]/(double)Nx;
    mom[i+2*Vol]=2.*mPI*x[2]/(double)Ny;
    
    latmom[i]=2.*sin(mom[i]/2.);
    latmom[i+Vol]=2.*sin(mom[i+Vol]/2.);
    latmom[i+2*Vol]=2.*sin(mom[i+2*Vol]/2.);
    
  }
  
  
  int t_fl,index;
  double sum;
  double scalarprod;
  double latmomsq;
  
    for(t_fl=0;t_fl<=N_steps;t_fl++){
      
      sum=0.;
      
      for(index=0;index<Vol;index++){
	scalarprod=mom[index]*x1+mom[index+Vol]*x2+mom[index+2*Vol]*x3;
	latmomsq=latmom[index]*latmom[index]+latmom[index+Vol]*latmom[index+Vol]+latmom[index+2*Vol]*latmom[index+2*Vol];
	sum+=-latmomsq*cos(scalarprod)*exp(-t_fl*dt*latmomsq);
      }
      
      sum/=(double)Vol;
      dt_jac_v[t_fl]=sum;
      
    }

    
  free(latmom);
  free(mom);
  
  mylog("[JACOBIAN]: Tabulation ended \n");
  mylog("\n");
  
    return;
}

    



  
  /* FREE_J: function used for freeing Jacobian vector memory */

  
  
  void free_J(double *J){

  free(J);
  J=NULL;

}





/* FREE_JAC: function for freeing the Jacobian table  */



void free_jac(void){

  free(jac_v);
  return;

}


void free_dt_jac(void){

  free(dt_jac_v);
  return;

}
