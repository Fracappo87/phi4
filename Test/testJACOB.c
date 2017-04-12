/* JACOBTEST: library used for testing Runge Kutta integrators for the Jacobian matrix */



#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include"testJACOB.h"
#include"Wilson_Flow.h"
#include"Jacobian.h"
#include"geometry.h"
#include"lattice_config.h"





/* DRHO_PHIANALYTICCOS: function used to define the analytic solution to be compared with the numerical one                                        */
/* In our case, the analytic initial condition is the rho-direction derivative of a real component of a plane wave of momentum "p" (cos(px))       */
/* The evolved solution can be trivially shown to be -(p_L_rho)sin(p(x+rho/2))exp[(-(p_L)^2+m^2)t] where p_L is the lattice momentum               */
/* -define the exponent of the heat kernel at finite lattice spacing and given bare mass of the flow equation                                      */
/* L-for every dimension i                                                                                                                         */
/*  -define the product xp                                                                                                                         */
/* -define the evolved solution as given above                                                                                                     */



double drho_phianalyticcos(double t,int i,int rho, double *p){
  
  int x[d];
  double xp=0.0;
  double exponent=-t*(RK.flow_mass*RK.flow_mass+sq_lattice_mom(p));
  
  coordinates(i,x);
  
  for(i=0;i<d;i++) {
   
    if(i==rho) xp+=p[i]*(x[i]+1/2);
    else xp+=p[i]*x[i];
  
  } 
 
  return -lattice_mom(p,rho)*exp(exponent)*sin(xp); 
  
}





/* TESTJACOB_INIT: function used to set the initial condition of the numerical solution. */



void testJACOB_init(double *p, double *X, double *J){

  int i;
  
  for(i=0;i<Vol;i++) X[i] = phianalyticcos(0, i, p);
  for(i=0;i<Vol;i++) J[i] = drho_phianalyticcos(0, i, 0, p);
  
}





/* TESTJACOB: function used to test the efficency of a given integration algorithm                                                                                                  */
/* The purpose of the test is measuring the difference between analytic and numerical solution                                                                                       */
/* If the test is correct, the difference should scale with a power law, with an exponent equal to the order of the integrator.                                                      */
/* As initial condition, we will use the rho-derivative of the monochromatic wave f_p(x)=cos(xp): the analytic solution can be esaily evaluated from the discretized flow equation   */
/* The test starts with some trivial checks on the volume mean of the squared initial condition which, if normalized by a factor 2/V, should be equal to 1                           */
/* Then a rapid check is made on the value of the lattice squared moment components                                                                                                  */
/* Starting the test                                                                                                                                                                 */
/* L-for each integration step                                                                                                                                                       */ 
/*  -set to zero the sum variables                                                                                                                                                   */
/*  L-for each lattice site                                                                                                                                                          */
/*   -define the volume mean of the squared analytic and numerical solution                                                                                                          */
/*   -define the L2 norm of their difference                                                                                                                                         */
/*  -renormalize each sum variable by the appropriate volume factor, print it in the logfile                                                                                         */
/*  -make an integration step                                                                                                                                                        */
/*  IF-the order of the integrator is equal to 1: make an Euler integration step                                                                                                     */
/*  ELSE IF-the order of the integrator is equal to 2: make a 2nd order RK integration step                                                                                          */
/*  ELSE IF-the order of the integrator is equal to 3: make a 3rd order RK integration step                                                                                          */
/*  ELSE IF-the order of the integrator is equal to 4: make a 4th order RK integration step                                                                                          */
/* -free vectors, finish the test                                                                                                                                                    */



void testJACOB(){

  int i,rho=2;
  double *Xi_s,*J_s;
  int t=0;
  double sum4=0.;  
  double sum5=0.0;
  double Norm2diff=0.;
  double p[d];
  double Vol_mean=0.0;
  double Vol_mean2=0.0;
  double Vol_meanJ=0.0;
  double Vol_mean2J=0.0;
  
  p[0]=2.*mPI/(double)Nt;
  p[1]=2.*mPI/(double)Nx;
  p[2]=2.*mPI/(double)Ny;
  
  Xi_s=(double*)malloc(Vol*sizeof(double));
  J_s=(double*)malloc(Vol*sizeof(double));
  Wilson_Flow_init(); 
  testJACOB_init(p,Xi_s,J_s);
  
  mylog("#################################################################\n");
  mylog("#                                                               #\n");
  mylog("#          PROGRAM FOR TESTING RUNGE KUTTA INTEGRATOR           #\n");
  mylog("#                OF THE JACOBIAN FLOW EQUATION                  #\n");
  mylog("#                                                               #\n");
  mylog("#################################################################\n");
  mylog("\n");

  for(i=0;i<Vol;i++){
    
    Vol_mean+=Xi_s[i];
    Vol_mean2+=Xi_s[i]*Xi_s[i];
    Vol_meanJ+=J_s[i];
    Vol_mean2J+=J_s[i]*J_s[i];
  }

  mylog("[TESTJACOB]: The purpose of the test is measuring the difference between analytic and numerical solution \n");
  mylog("[TESTJACOB]: If the test is correct, the difference should scale with a power law, with an exponent equal to the order of the integrator \n");
  mylog("[TESTJACOB]: We will test the RK%d integrator \n",int_ord);
  mylog("[TESTJACOB]: As initial condition, we will use the monochromatic wave f_p(x)=cos(xp): the analytic solution can be easily evaluated from the discretized flow equation \n");
  mylog("[TESTJACOB]: It should be stressed that the solution given by our Runge Kutta integrator represent the action of the Jacobian matrix on a flowed field at a given time step t* \n");
  mylog("[TESTJACOB]: Keeping this in mind, we must remember that the analytic solution chosen has to represent the evolution of the action of such Jacobian, with initial condition f_p(x)=cos(xp) \n");
  mylog("[TESTJACOB]: Starting with some trivial checks on the initial condition \n");
  mylog("\n");
  mylog("[TESTJACOB]: Time p-component p[0] = %lf \n",p[0]);
  mylog("[TESTJACOB]: Space p-component p[1] = %lf \n",p[1]);
  mylog("[TESTJACOB]: Space p-component p[2] = %lf \n",p[2]);
  mylog("\n");
  mylog("[TESTJACOB]: Starting with the measurement \n");
  mylog("\n");
  mylog("[TESTJACOB]: Volume mean initial field condition = %e \n",Vol_mean*(2/(double)Vol));
  mylog("[TESTJACOB]: Volume mean^2 initial field condition = %e \n",Vol_mean2*(2/(double)Vol));
  mylog("[TESTJACOB]: Volume mean initial Jacobian-field condition = %e \n",Vol_meanJ*(2/(double)Vol));
  mylog("[TESTJACOB]: Volume mean^2 initial Jacobian-field condition = %e \n",Vol_mean2J*(2/(double)Vol));
  
 for(i=0;i<d;i++){

  mylog("[TESTJACOB]: Lattice mom check = %lf \n",lattice_mom(p,i));

  }

  mylog("[TESTJACOB]: Square lattice mom check = %lf \n",sq_lattice_mom(p));
 mylog("[TESTJACOB]: We check the solution for the %d component \n",rho);
  mylog("\n");
  
  for(t=0;t<N_steps+1;t++){
    
    sum4=0.0;
    sum5=0.0;    
    Norm2diff=0.0;
    
    for(i=0;i<Vol;i++){
      
      sum4+=J_s[i]*J_s[i];
      sum5+=drho_phianalyticcos(t*dt, i, rho, p)*drho_phianalyticcos(t*dt, i, rho, p);
      Norm2diff+=(J_s[i]-drho_phianalyticcos(t*dt, i, rho, p))*(J_s[i]-drho_phianalyticcos(t*dt, i, rho, p));
      
    }
    
    mylog("[TESTJACOB]: Time = %f -- Exact solution = %e \n",t*dt,sum5*(2.0)/(double)Vol);
    mylog("[TESTJACOB]: Time = %f -- Num solution = %e  \n",t*dt,sum4*(2.0)/(double)Vol);
    mylog("[TESTJACOB]: Time = %f -- Diff = %e  \n",t*dt,sqrt(Norm2diff)/(double)Vol);
    mylog("\n");

    if(int_ord==1) Jacobian_Flow_RK1(dt,Xi_s,J_s);
    else if (int_ord==2) Jacobian_Flow_RK2(dt,Xi_s,J_s);
    else if (int_ord==3) Jacobian_Flow_RK3(dt,Xi_s,J_s);
    else if (int_ord==4) Jacobian_Flow_RK4(dt,Xi_s,J_s);

    if(int_ord==1) Wilson_Flow_RK1(dt,Xi_s);
    else if (int_ord==2) Wilson_Flow_RK2(dt,Xi_s);
    else if (int_ord==3) Wilson_Flow_RK3(dt,Xi_s);
    else if (int_ord==4) Wilson_Flow_RK4(dt,Xi_s);
    
  }
  
  free_WF(Xi_s);
  free_J(J_s);  

}


