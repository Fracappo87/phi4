/*MAIN PROGRAM: it does what I want */



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





/* BEGINNING */



int main(int argc, char *argv[]){
  
  int i,j,k,t,count;
  int counter=0;
  int ncnfg;
  double DeltaT;
  double accettanza=0., accettanza_cluster=0.;
  struct timeval t1,t2;
  struct timeval T1,T2;
  //  int invdt;
  
  
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
  
  
  
  
  
  /* READING FLOW TIMES */  



  flow_input("flowtimes.dat");





  /* INITIALIZING RANDOM NUMBER GENERATOR */
  
  
  
  rlxd_init(2,seed);
  
  
  
  
  
  /* GEOMETRY SETTING */
  
  
  
  geometry();
  
  
  
  
  
  /* SETTING INITIAL FIELDS CONFIGURATION */
  
  
  
  configuration();
  
  
  
  
  
  /* SETTING INITAL CLUSTER CONFIGURATION */
  
  
  
  cluster_init();
  

  mylog("[MAIN]: Observables measured at the following flow times \n");
  
  for(k=0;k<N_flowtimes;k++){
    
    mylog("[MAIN]: t = %e ; c(t)= %e\n",Flow_vec[k]*dt,sqrt(Flow_vec[k]*dt*8)/((double)Nt));
    
  }
  
  mylog("\n");
  
  
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
  
  

  
  
  /* MEASUREMENT */
  
  
  
  mylog("[MAIN]: Starting with the measurements\n");

  double *Xi=(double*)malloc(Vol*sizeof(double));

  Wilson_Flow_init();

  for(k=0;k<N_traj;k++){
    
    counter++;
    gettimeofday(&t1, NULL);
    
    for(i=0;i<N_meas;i++){
      
      update_sweep(N_sweeps,sweep_ratio);

      mylog("[MAIN]: Couplings differing order a^2 \n");
      mylog("[MAIN]: <phi4> %lf %e \n",0.,phi4(Phi));
      mylog("[MAIN]: <phi4ps> %lf %e \n",0.,corr(1,Phi));	  
      
	  
      
      
      for(j=0;j<Vol;j++) {
	
	Xi[j]=Phi[j];
	
      }
      
      

      count=0;
      
      for(t=0;t<N_steps;t++){
	
	Wilson_Flow_RK4(dt,Xi);
	
	if((Flow_vec[count]-t-1)==0){
	      
	  mylog("[MAIN]: Couplings differing order a^2 \n");
	  mylog("[MAIN]: <phi4> %lf %e \n",(t+1)*dt,phi4(Xi));
	  mylog("[MAIN]: <phi4ps> %lf %e \n",(t+1)*dt,corr(1,Xi));	  
	    count ++;
	    
	}
      }   
    }
    
    
    gettimeofday(&t2, NULL);
    
    DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
    mylog("[MAIN]: Observables measured in %lf seconds\n", accettanza);
    
    /* if(counter%N_report==0){
       
       DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
       accettanza=(double)acc/call;
       accettanza_cluster=(double)acc_clu/call_clu;
       ncnfg=(k+1)+N_term+Last_count;
       
       mylog("[MAIN] Configuration %d generated in %lf seconds\n", ncnfg,(double)DeltaT/1000000.0);
       mylog("[MAIN] Metropolis acceptance = %lf\n", accettanza);
       mylog("[MAIN] Cluster acceptance =  %lf\n", accettanza_cluster);
       mylog("\n");
       
       } */    
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
  free(Flow_vec);  
  mylog("[MAIN]: closing program\n");
  
  return 0;
  
}
