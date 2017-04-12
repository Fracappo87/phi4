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
#define INTEGRATED
#include"TWI_obs_engine.h"




/* BEGINNING */



int main(int argc, char *argv[]){
  
  int i,j,k,t,count;
  int counter=0;
  int ncnfg;
  double DeltaT;
  double accettanza=0., accettanza_cluster=0.;
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
  
  
  
  
  
  /* READING FLOW TIMES */  



  flow_input("flowtimes.dat");





  /* READING PROBE POSITIONS */



  TWI_input("twi.dat");





  /* INITIALIZING RANDOM NUMBER GENERATOR */
  
  
  
  rlxd_init(2,seed);
  
  
  
  
  
  /* GEOMETRY SETTING */
  
  
  
  geometry();
  
  
  
  
  
  /* SETTING INITIAL FIELDS CONFIGURATION */
  
  
  
  configuration();
  
  
  
  
  
  /* SETTING INITAL CLUSTER CONFIGURATION */
  
  
  
  cluster_init();







  /*  FLOW TIMES INFO  */

  mylog("[MAIN]: Observables measured at the following flow times \n");
  
  for(k=0;k<N_flowtimes;k++){
    
    mylog("[MAIN]: t = %e ; c(t)= %e\n",Flow_vec[k]*dt,sqrt(Flow_vec[k]*dt*8)/((double)Nt));
    
  }
  
  mylog("\n");





  /* TWI INFO */


  if(lhs_measure)  
  mylog("[MAIN]: Measuring the LHS of the TWI \n");
  if(rhs_measure)  
  mylog("[MAIN]: Measuring the RHS of the TWI \n");
  mylog("\n");


  double *avg_vector=avg_vector_init();

  vec_clearer(avg_vector);


  mylog("[MAIN]: Measurements averaged over intermediate blocks of size %d \n",avg_TWI);
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

  gettimeofday(&t1, NULL);
  jac();
  gettimeofday(&t2, NULL);

  DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
  mylog("[MAIN]: Jacobian tabulation done in %lf seconds \n",(double)DeltaT/1000000.0);
  mylog("\n");



  for(k=0;k<N_traj;k++){
    
    counter++;
    gettimeofday(&t1, NULL);
    
    for(i=0;i<N_meas;i++){
      
      update_sweep(N_sweeps,sweep_ratio);

      int np,ntop;


      
      if(lhs_measure){
	for(np=0;np<NUM_PROB;np++){
	  //	  mylog("[MAIN]: <delta_V[%d]_test> %lf %e \n",np+1,0.,LHS_delta_V(np+1,traslpos,probepos,0,Phi,Phi));    
	  avg_vector[np]+=LHS_delta_V_int(np+1,0,Phi,Phi);
	}
      }
      
      if(rhs_measure){	
	for(np=0;np<NUM_PROB;np++)
	  for(ntop=0;ntop<NUM_OP;ntop++){
	    //mylog("[MAIN]: <RHS[%d,%d]_test> %lf %e \n",np+1,ntop+1,0.,RHS(traslpos,probepos,np,ntop,Phi,Phi));    
	    avg_vector[NUM_PROB+np+NUM_PROB*ntop]+=RHS_int(np,ntop,Phi,Phi);
	    }
      }
      
      for(j=0;j<Vol;j++) {
	
	Xi[j]=Phi[j];
	
      }
      
      

      count=0;
      
      for(t=0;t<N_steps;t++){
	
	Wilson_Flow_RK4(dt,Xi);
	
	if((Flow_vec[count]-t-1)==0){

	  if(lhs_measure){	  
	    for(np=0;np<NUM_PROB;np++){
	      //     mylog("[MAIN]: <delta_V[%d]_test> %lf %e \n",np+1,(t+1)*dt,LHS_delta_V(np+1,traslpos,probepos,t+1,Xi,Phi));
	      avg_vector[NUM_PROB*(NUM_OP+1)*(count+1)+np]+=LHS_delta_V_int(np+1,t+1,Xi,Phi);
	    }
	  }
	  
	  if(rhs_measure){
	    
	    for(np=0;np<NUM_PROB;np++)
	      for(ntop=0;ntop<NUM_OP;ntop++){
		//mylog("[MAIN]: <RHS[%d,%d]_test> %lf %e \n",np+1,ntop+1,(t+1)*dt,RHS(traslpos,probepos,np,ntop,Xi,Phi));    
		avg_vector[NUM_PROB+np+NUM_PROB*ntop+NUM_PROB*(NUM_OP+1)*(count+1)]+=RHS_int(np,ntop,Xi,Phi);
	  }	    
	}
	  count ++;
	  
	}
      }
      
      
      if((i+1)%avg_TWI==0){

      mylog("[MAIN]: Probe observables variation \n");
	
	if(lhs_measure){
	  for(np=0;np<NUM_PROB;np++)
	    mylog("[MAIN]: <delta_V[%d]> %lf %e \n",np+1,0.,avg_vector[np]/(double)avg_TWI);                                                                                                         
	    
	}

	if(rhs_measure){
	  for(np=0;np<NUM_PROB;np++)
	    for(ntop=0;ntop<NUM_OP;ntop++)

	      mylog("[MAIN]: <RHS[%d,%d]> %lf %e \n",np+1,ntop+1,0., avg_vector[NUM_PROB+np+NUM_PROB*ntop]/(double)avg_TWI);                                                                                                       	     
	}

	
	for(t=0;t<N_flowtimes;t++){
	  
	  if(lhs_measure){
            for(np=0;np<NUM_PROB;np++)
	      
	      mylog("[MAIN]: <delta_V[%d]> %lf %e \n",np+1,Flow_vec[t]*dt,avg_vector[NUM_PROB*(NUM_OP+1)*(t+1)+np]/(double)avg_TWI);                                                                                            	    
	  }
	  
	  if(rhs_measure){
	    
	    for(np=0;np<NUM_PROB;np++)
	      for(ntop=0;ntop<NUM_OP;ntop++)
	
		mylog("[MAIN]: <RHS[%d,%d]> %lf %e \n",np+1,ntop+1,Flow_vec[t]*dt,avg_vector[NUM_PROB+np+NUM_PROB*ntop+NUM_PROB*(NUM_OP+1)*(t+1)]/(double)avg_TWI);                                                                                               
	  }
	}

	vec_clearer(avg_vector);      
	
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
  free(avg_vector);
  mylog("[MAIN]: closing program\n");
  
  return 0;
  
}
