/* MAIN_TWI: main program for the evaluation of the translation operator renormalization constant */



#define MAIN_PROGRAM

#include<stdio.h>
#include <time.h>
#include <sys/time.h>

#include"geometry.h"
#include"glob_const.h"
#include"input_output.h"
#include"lattice_config.h"
#include"update.h"
#include"observables.h"
#include"Wilson_Flow.h"
#include"Jacobian.h"




/* BEGINNING */



int main(int argc, char *argv[]){

  int i,j,k,t;
  int counter=0;
  int ncnfg;
  double DeltaT;
  double accettanza, accettanza_cluster;
  struct timeval t1,t2;
  struct timeval T1,T2;





  /* PRESENTATION */
  

  
  init_mylog(argv[1],argc);
  
  gettimeofday(&T1, NULL); 
  
  mylog("####################################################\n");
  mylog("#                                                  #\n");
  mylog("#          PROGRAM FOR SCALAR PHI THEORY           #\n");
  mylog("#                                                  #\n");
  mylog("####################################################\n");
  mylog("\n");
  
  
  
  
  
  /* READING INPUT PARAMETERS*/


 
  readinput("parameters.dat");





  /* PROGRAM'S PARAMETERS CHECK */
  
  
  
  mylog("[MAIN]: Number of trajectories = %d\n",N_traj);
  mylog("[MAIN]: Number of sweeps for each measure, N_sweeps = %d\n",N_sweeps);
  mylog("[MAIN]: Metropolis/Cluster sweep ratio, sweep_ratio = %d\n",sweep_ratio);
  mylog("[MAIN]: Number of thermalization steps, N_term = %d\n",N_term);
  mylog("[MAIN]: Report rate, N_report = %d\n",N_report);
  mylog("[MAIN]: Save rate, N_save = %d\n",N_save);
  mylog("[MAIN]: Measure rate, N_save = %d\n",N_meas);
  mylog("\n");
  mylog("[MAIN]: UPDATE ALGORITHM PARAMETERS\n");
  mylog("[MAIN]: Metropolis delta parameter = %lf\n",delta);
  mylog("[MAIN]: Random seed = %d\n",seed);
  mylog("[MAIN]: Mass parameter, mass = %lf\n",mass);
  mylog("[MAIN]: Quartic interaction coupling, lambda = %lf\n",lambda);
  mylog("[MAIN]: Cubic interaction coupling, g = %lf\n",g);
  mylog("[MAIN]: Sixth interaction coupling, kappa = %lf\n",kappa);

  if(mass_square>0){
    
    mylog("[MAIN]: Positive square mass simulation\n");   
    mylog("[MAIN]: Square mass = %lf\n",mass_square);
    mylog("\n");

  }else{
    
    mylog("[MAIN]: Negative square mass simulation\n");
    mylog("[MAIN]: Square mass = %lf\n",mass_square);
    mylog("\n");
    
  } 
  
  
  
  
  
  /* INITIALIZING RANDOM NUMBER GENERATOR */


 
  rlxd_init(2,seed);
  
  
  
  
  
  /* GEOMETRY SETTING */
  
  
  
  geometry();
  
  
  
  
  
  /* SETTING INITIAL FIELDS CONFIGURATION */
  


  configuration();




  
  /* SETTING INITAL CLUSTER CONFIGURATION */
  
  
  
  cluster_init();


  
 
 
  /* THERMALIZATION */
  
  
  
  mylog("[MAIN] Starting thermalization\n");
  mylog("\n");
  
  for(k=0;k<N_term/N_report;k++){
    for(i=0;i<N_report-1;i++){
      
      update_sweep(N_sweeps,sweep_ratio); 
      
    }
    
    gettimeofday(&t1, NULL);
    update_sweep(N_sweeps,sweep_ratio); 
    gettimeofday(&t2, NULL);
    
    DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
    accettanza=(double)acc/call;
    accettanza_cluster=(double)acc_clu/call_clu;
    ncnfg=(k+1)*N_report+Last_count;
    
    mylog("[MAIN] Configuration %d generated in %lf seconds\n", ncnfg,(double)DeltaT/1000000.0);
    mylog("[MAIN] Metropolis acceptance = %lf\n", accettanza);
    mylog("[MAIN] Cluster acceptance =  %lf\n", accettanza_cluster);
    mylog("\n");
    
  }
   
  
 
   
  
  /* MEASUREMENTS  */
  /* Zdelta measurement: according to the prescription given for the evaluation of Z_delta, we need to take two flowe probe observables at two different points along the time (y0,x0) direction    */
  /* The first of two (y0) must lie outside the integration domani centered in (x0),[x0-l,x0+l]                                                                                                     */
  /* For every configuration simulated, we set the initial conditon for the field and thsssse integral action of the Jacobian operator, then evolve for N_steps Wilson Flow steps                   */
  /* It should be underlined that the intial condition for the Jacobian needs an ad hoc customization, in order to effectively reduce its integral action to the [x0-l,x0+l]x[0,L]^{d-1} subdomain  */



  mylog("[MAIN]: Starting with the measurements\n");
  
 
  //double *Jf=(double*)malloc(Vol*sizeof(double));
  double *Xi=(double*)malloc(Vol*sizeof(double));
  double *dPhi=(double*)malloc(Vol*sizeof(double));
 
  
  
  int xext=0;
  int yint=Nt/2;
  int l=Nt/4;  
  
  mylog("[MAIN]:Integration slice : [%d-%d,%d+%d] \n", yint,l,yint,l);      
  mylog("[MAIN]:External probe point : [%d] \n", xext);
  

  Wilson_Flow_init();

  gettimeofday(&t1, NULL);
  jac(0,0,0);
  gettimeofday(&t2, NULL);

  DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
  mylog("[MAIN]: Jacobian tabulation done in %lf seconds \n",(double)DeltaT/1000000.0);
  mylog("\n");



  for(k=0;k<N_traj;k++){
    
    counter++;
    gettimeofday(&t1, NULL);
    
    for(i=0;i<N_meas;i++){
      
      update_sweep(N_sweeps,sweep_ratio);


	
      
      for(j=0;j<Vol;j++) {
	
	Xi[j]=Phi[j];
	//Jf[j]=(Phi[d_up[j*d]]-Phi[j]);
	dPhi[j]=(Phi[d_up[j*d]]-Phi[j]);
	
      }

      
      //Zdelta_set_condition(yint,l,Jf);
      

      mylog("[MAIN]:Integrated T.W.I: (flowtime/Observable) \n");
      //      mylog("[MAIN]: <L.H.S.f> %lf %e \n",0.,LHpert(xext,yint,Phi,Jf));      	
      mylog("[MAIN]: <L.H.S.f.exact> %lf %e \n",0.,LHpert_exact(xext,yint,l,0,dPhi,Xi,jac_v));      	
      mylog("[MAIN]: <L.H.S.f.exact.vol> %lf %e \n",0.,LHpert_exact_vol(yint-xext,l,0,dPhi,Xi,jac_v));      	
      mylog("[MAIN]: <R.H.S.f.naive> %lf %e \n",0.,corr(yint-xext+1,Xi)-corr(yint-xext,Xi));      	
      mylog("[MAIN]: <R.H.S.f.chain> %lf %e \n",0.,-rhschainnew(yint-xext,Xi));      	
	

      
      for(t=0;t<N_steps;t++){
	
	//Jacobian_Flow_RK4(dt,Xi,Jf);
	Wilson_Flow_RK4(dt,Xi);
	
	
	if((t+1)%N_flow==0){
	  
	  mylog("[MAIN]:Integrated T.W.I: (flowtime/Observable) \n");
	  //mylog("[MAIN]: <L.H.S.f> %lf %e \n",(1+t)*dt,LHpert(xext,yint,Xi,Jf));  
	  mylog("[MAIN]: <L.H.S.f.exact> %lf %e \n",(1+t)*dt,LHpert_exact(xext,yint,l,1+t,dPhi,Xi,jac_v));      	
	  mylog("[MAIN]: <L.H.S.f.exact.vol> %lf %e \n",(1+t)*dt,LHpert_exact_vol(yint-xext,l,1+t,dPhi,Xi,jac_v));      	      	    	    	
	  mylog("[MAIN]: <R.H.S.f.naive> %lf %e \n",(1+t)*dt,corr(yint-xext+1,Xi)-corr(yint-xext,Xi));      	
	  mylog("[MAIN]: <R.H.S.f.chain> %lf %e \n",(1+t)*dt,-rhschainnew(yint-xext,Xi));      	
	}
	
      }
      
      gettimeofday(&t2, NULL);
      
      DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
      
      if(counter%N_report==0){
	
	DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
	accettanza=(double)acc/call;
	accettanza_cluster=(double)acc_clu/call_clu;
	ncnfg=(k+1)+N_term+Last_count;
	
	mylog("[MAIN] Configuration %d generated in %lf seconds\n", ncnfg,(double)DeltaT/1000000.0);
	mylog("[MAIN] Metropolis acceptance = %lf\n", accettanza);
	mylog("[MAIN] Cluster acceptance =  %lf\n", accettanza_cluster);
	mylog("\n");
	
      }     
    }
  }  
  
  gettimeofday(&T2, NULL);
  
  DeltaT=(T2.tv_sec*1000000+T2.tv_usec)-(T1.tv_sec*1000000+T1.tv_usec);
  mylog("[MAIN]: Total simulation time = %lf seconds \n",(double)DeltaT/1000000.0);
  
  
    
  
  
  /* FREEING MEMORY/CLOSING OUTPUT FILES*/ 
  
    
  
  mylog("[MAIN]: Freeing memory, closing files");
  mylog("\n");
  
  free_geo();
  free_cluster();
  free_Lattice(); 
  free_WF(Xi);
  //free_J(Jf);
  free_jac();
  
  mylog("[MAIN]: closing program\n");
  
  return 0;
  
}  


