/* RKTEST: library used for testing Runge Kutta integrators */



#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include"testRK.h"
#include"Wilson_Flow.h"
#include"geometry.h"
#include"lattice_config.h"





/* LATTICE_MOM: trivial function we need in order to define the lattice momentum (see any lattice field theory reference for more infos )*/



double lattice_mom(double *p,int mu){

  return 2*sin(p[mu]/2);

}





/* SQ_LATTICE_MOM: trivial function we need in order to define the square lattice momentum (see any lattice field theory reference for more infos )*/



double sq_lattice_mom(double *p){
  
  int i;
  double sq_lat_mom=0;
  
  for(i=0;i<d;i++) sq_lat_mom+=4*sin(p[i]/2)*sin(p[i]/2);
  
  return sq_lat_mom;  

}





/* PHIANALYTICCOS: function used to define the analytic solution to be compared with the numerical one               */
/* In our case, the analytic initial condition is the real component of a plane wave of momentum "p" (cos(px))       */
/* The evolved solution can be trivially shown to be cosp(px)exp[(-(p_L)^2+m^2)t] where p_L is the lattice momentum  */
/* -define the exponent of the heat kernel at finite lattice spacing and given bare mass of the flow equation        */
/* L-for every dimension i                                                                                           */
/*  -define the product xp                                                                                           */
/* -define the evolved solution as given above                                                                       */



double phianalyticcos(double t,int i, double *p){
  
  int x[d];
  double xp=0.0;
  double exponent=-t*(RK.flow_mass*RK.flow_mass+sq_lattice_mom(p));
  
  coordinates(i,x);
  
  for(i=0;i<d;i++) xp+=p[i]*x[i];
  
  return exp(exponent)*cos(xp);
  
}





/* TESTRK_INIT: function used to set the initial condition of the numerical solution. */



void testRK_init(double *p, double *X){

  int i;
  
  for(i=0;i<Vol;i++) X[i] = phianalyticcos(0, i, p);
  
}





/* TESTRK: function used to test the efficency of a given integration algorithm                                                                                */
/* The purpose of the test is measuring the difference between analytic and numerical solution                                                                 */
/* If the test is correct, the difference should scale with a power law, with an exponent equal to the order of the integrator.                                */
/* As initial condition, we will use the monochromatic wave f_p(x)=cos(xp): the analytic solution can be esaily evaluated from the discretized flow equation   */
/* The test starts with some trivial checks on the volume mean of the squared initial condition which, if normalized by a factor 2/V, should be equal to 1     */
/* Then a rapid check is made on the value of the lattice squared moment components                                                                            */
/* Starting the test                                                                                                                                           */
/* L-for each integration step                                                                                                                                 */ 
/*  -set to zero the sum variables                                                                                                                             */
/*  L-for each lattice site                                                                                                                                    */
/*   -define the volume mean of the squared analytic and numerical solution                                                                                    */
/*   -define the L2 norm of their difference                                                                                                                   */
/*  -renormalize each sum variable by the appropriate volume factor, print it in the logfile                                                                   */
/*  -make an integration step                                                                                                                                  */
/*  IF-the order of the integrator is equal to 1: make an Euler integration step                                                                               */
/*  ELSE IF-the order of the integrator is equal to 2: make a 2nd order RK integration step                                                                    */
/*  ELSE IF-the order of the integrator is equal to 3: make a 3rd order RK integration step                                                                    */
/*  ELSE IF-the order of the integrator is equal to 4: make a 4th order RK integration step                                                                    */
/* -free vectors, finish the test                                                                                                                              */



void testRK(){

  int i;
  double *Xi_s;
  int t=0;
  double sum4=0.;  
  double sum5=0.0;
  double Norm2diff=0.;
  double p[d];
  double Vol_mean=0.0;
  double Vol_mean2=0.0;
  
  p[0]=2.*mPI/(double)Nt;
  p[1]=2.*mPI/(double)Nx;
  p[2]=2.*mPI/(double)Ny;
  
  Xi_s=(double*)malloc(Vol*sizeof(double));
  Wilson_Flow_init(); 
  testRK_init(p,Xi_s);
  
  mylog("#################################################################\n");
  mylog("#                                                               #\n");
  mylog("#          PROGRAM FOR TESTING RUNGE KUTTA INTEGRATOR           #\n");
  mylog("#                 OF THE FIELD FLOW EQUATION                    #\n");
  mylog("#                                                               #\n");
  mylog("#################################################################\n");
  mylog("\n");

  for(i=0;i<Vol;i++){
    
    Vol_mean+=Xi_s[i];
    Vol_mean2+=Xi_s[i]*Xi_s[i];
    
  }

  mylog("[TESTRK]: The purpose of the test is measuring the difference between analytic and numerical solution \n");
  mylog("[TESTRK]: If the test is correct, the difference should scale with a power law, with an exponent equal to the order of the integrator \n");
  mylog("[TESTRK]: We will test the RK%d integrator \n",int_ord);
  mylog("[TESTRK]: As initial condition, we will use the monochromatic wave f_p(x)=cos(xp): the analytic solution can be easily evaluated from the discretized flow equation \n");
  mylog("[TESTRK]: Starting with some trivial checks on the initial condition \n");
  mylog("\n");
  mylog("[TESTRK]: Time p-component p[0] = %lf \n",p[0]);
  mylog("[TESTRK]: Space p-component p[1] = %lf \n",p[1]);
  mylog("[TESTRK]: Space p-component p[2] = %lf \n",p[2]);
  mylog("\n");
  mylog("[TESTRK]: Starting with the measurement \n");
  mylog("\n");
   mylog("[TESTRK]: Volume mean initial condition = %e \n",Vol_mean*(2/(double)Vol));
  mylog("[TESTRK]: Volume mean^2 initial condition = %e \n",Vol_mean2*(2/(double)Vol));

  for(i=0;i<d;i++){

  mylog("[TESTJACOB]: Lattice mom check = %lf \n",lattice_mom(p,i));

  }
  
  mylog("[TESTRK]: Square lattice mom check = %lf \n",sq_lattice_mom(p));
  mylog("\n");
  
  for(t=0;t<N_steps+1;t++){
    
    sum4=0.0;
    sum5=0.0;    
    Norm2diff=0.0;
    
    for(i=0;i<Vol;i++){
      
      sum4+=Xi_s[i]*Xi_s[i];
      sum5+=phianalyticcos(t*dt, i, p)*phianalyticcos(t*dt, i, p);
      Norm2diff+=(Xi_s[i]-phianalyticcos(t*dt, i, p))*(Xi_s[i]-phianalyticcos(t*dt, i, p));
      
    }
    
    mylog("[TESTRK]: Time = %f -- Exact solution = %e \n",t*dt,sum5*(2.0)/(double)Vol);
    mylog("[TESTRK]: Time = %f -- Num solution = %e  \n",t*dt,sum4*(2.0)/(double)Vol);
    mylog("[TESTRK]: Time = %f -- Diff = %e  \n",t*dt,sqrt(Norm2diff)/(double)Vol);
    mylog("\n");

    if(int_ord==1) Wilson_Flow_RK1(dt,Xi_s);
    else if (int_ord==2) Wilson_Flow_RK2(dt,Xi_s);
    else if (int_ord==3) Wilson_Flow_RK3(dt,Xi_s);
    else if (int_ord==4) Wilson_Flow_RK4(dt,Xi_s);
    
  }
  
  free_WF(Xi_s);
  
}







