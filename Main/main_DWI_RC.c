/*MAIN PROGRAM: it does what I want */



#define MAIN_PROGRAM
//#define INTEGRATED
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
#include"DWI_obs_engine.h"




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
  /* Translation application point */
  /* Probe position                */



  int traslcoor[d];
  int probecoor[d];

  coordinates(traslpos,traslcoor);
  coordinates(probepos,probecoor);

  //  mylog("[MAIN]: Translation applied at the following point \n");
  // mylog("[MAIN]: index x=%d: vecx=[%d,%d,%d] \n",traslpos,traslcoor[0],traslcoor[1],traslcoor[2]);
  mylog("[MAIN]: Probe positioned at the following point \n");
  mylog("[MAIN]: index y=%d: vecy=[%d,%d,%d]\n",probepos,probecoor[0],probecoor[1],probecoor[2]);
  if(lhs_measure)  
  mylog("[MAIN]: Measuring the LHS of the DWI \n");
  if(rhs_measure)  
  mylog("[MAIN]: Measuring the RHS of the DWI \n");
  mylog("\n");


  double *avg_vector=dwiavg_vector_init();

  dwivec_clearer(avg_vector);


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
  
#ifdef INTEGRATED
  int volloop;
#endif  
  
  mylog("[MAIN]: Starting with the measurements\n");

  double *Xi=(double*)malloc(Vol*sizeof(double));

  Wilson_Flow_init();

  gettimeofday(&t1, NULL);
  dt_jac();
  gettimeofday(&t2, NULL);

  DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
  mylog("[MAIN]: Jacobian flow derivative, tabulation done in %lf seconds \n",(double)DeltaT/1000000.0);
  mylog("\n");


  //testing jacobian derivative
  /*  int vec[d];

  for(k=0;k<Vol;k++){
   
    coordinates(k,vec);
    mylog(" Jacobian[%d,%d,%d;%lf] =%lf\n",vec[0],vec[1],vec[2],Flow_vec[3]*dt,dt_jac_v[k+Vol*Flow_vec[3]]);

    }*/

  for(k=0;k<N_traj;k++){
    
    counter++;
    gettimeofday(&t1, NULL);
    
    for(i=0;i<N_meas;i++){
      
      update_sweep(N_sweeps,sweep_ratio);

      int np,ntop;


      
      if(lhs_measure){
	for(np=0;np<dwiNUM_PROB;np++){
	  //	  mylog("[MAIN]: <delta_V[%d]_test> %lf %e \n",np+1,0.,LHS_delta_V(np+1,traslpos,probepos,0,Phi,Phi));    
#ifndef INTEGRATED
	  avg_vector[np]+=dwiLHS(np,probepos,0,Phi,Phi);
#else
	  
	  for(volloop=0;volloop<Vol;volloop++)
	    avg_vector[np]+=dwiLHS(np,volloop,0,Phi,Phi);
#endif
	}
      }
      
      if(rhs_measure){	
	for(np=0;np<dwiNUM_PROB;np++)
	  for(ntop=0;ntop<dwiNUM_OP;ntop++){
	    //mylog("[MAIN]: <RHS[%d,%d]_test> %lf %e \n",np+1,ntop+1,0.,RHS(traslpos,probepos,np,ntop,Phi,Phi));    
#ifndef INTEGRATED
	    avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop]+=dwiRHS(np,ntop,probepos,Phi,Phi);
#else
	    for(volloop=0;volloop<Vol;volloop++)
	      avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop]+=dwiRHS(np,ntop,volloop,Phi,Phi);
#endif

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
	    for(np=0;np<dwiNUM_PROB;np++){
	      //	      mylog("[MAIN]: <delta_V[%d]_test> %lf %e \n",np+1,(t+1)*dt,dwiLHS(np,probepos,t+1,Xi,Phi,Flow_vec));
#ifndef INTEGRATED
	      avg_vector[dwiNUM_PROB*(dwiNUM_OP+1)*(count+1)+np]+=dwiLHS(np,probepos,t+1,Xi,Phi);
#else
	      for(volloop=0;volloop<Vol;volloop++)
		avg_vector[dwiNUM_PROB*(dwiNUM_OP+1)*(count+1)+np]+=dwiLHS(np,volloop,t+1,Xi,Phi);
#endif

	    }
	  }
	  
	  if(rhs_measure){
	    
	    for(np=0;np<dwiNUM_PROB;np++)
	      for(ntop=0;ntop<dwiNUM_OP;ntop++){
		//mylog("[MAIN]: <RHS[%d,%d]_test> %lf %e \n",np+1,ntop+1,(t+1)*dt,RHS(traslpos,probepos,np,ntop,Xi,Phi));    
#ifndef INTEGRATED
		avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop+dwiNUM_PROB*(dwiNUM_OP+1)*(count+1)]+=dwiRHS(np,ntop,probepos,Xi,Phi);
#else
		for(volloop=0;volloop<Vol;volloop++)
		  avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop+dwiNUM_PROB*(dwiNUM_OP+1)*(count+1)]+=dwiRHS(np,ntop,volloop,Xi,Phi);
#endif
	      }	    
	}
	  count ++;
	  
	}
      }
      
      
      if((i+1)%avg_TWI==0){

      mylog("[MAIN]: Probe observables variation \n");
	
	if(lhs_measure){
	  for(np=0;np<dwiNUM_PROB;np++)
#ifndef INTEGRATED
	    mylog("[MAIN]: <delta_V[%d]> %lf %e \n",np+1,0.,avg_vector[np]/(double)avg_TWI);
#else
	  mylog("[MAIN]: <delta_V[%d]> %lf %e \n",np+1,0.,avg_vector[np]/(double)(avg_TWI*Vol));
#endif             
	    
	}

	if(rhs_measure){
	  for(np=0;np<dwiNUM_PROB;np++)
	    for(ntop=0;ntop<dwiNUM_OP;ntop++)
#ifndef INTEGRATED
	      mylog("[MAIN]: <RHS[%d,%d]> %lf %e \n",np+1,ntop+1,0., avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop]/(double)avg_TWI);                                                               
#else
	  mylog("[MAIN]: <RHS[%d,%d]> %lf %e \n",np+1,ntop+1,0., avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop]/(double)(avg_TWI*Vol));                                                           
#endif    
	}

	
	for(t=0;t<N_flowtimes;t++){
	  
	  if(lhs_measure){
            for(np=0;np<dwiNUM_PROB;np++)
	      
#ifndef INTEGRATED
	      mylog("[MAIN]: <delta_V[%d]> %lf %e \n",np+1,Flow_vec[t]*dt,avg_vector[dwiNUM_PROB*(dwiNUM_OP+1)*(t+1)+np]/(double)avg_TWI);                                                      
#else

	    mylog("[MAIN]: <delta_V[%d]> %lf %e \n",np+1,Flow_vec[t]*dt,avg_vector[dwiNUM_PROB*(dwiNUM_OP+1)*(t+1)+np]/(double)(avg_TWI*Vol));                                                      
#endif                                      	    
	  }
	  
	  if(rhs_measure){
	    
	    for(np=0;np<dwiNUM_PROB;np++)
	      for(ntop=0;ntop<dwiNUM_OP;ntop++)
#ifndef INTEGRATED	
		mylog("[MAIN]: <RHS[%d,%d]> %lf %e \n",np+1,ntop+1,Flow_vec[t]*dt,avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop+dwiNUM_PROB*(dwiNUM_OP+1)*(t+1)]/(double)avg_TWI);                  
#else
	    mylog("[MAIN]: <RHS[%d,%d]> %lf %e \n",np+1,ntop+1,Flow_vec[t]*dt,avg_vector[dwiNUM_PROB+np+dwiNUM_PROB*ntop+dwiNUM_PROB*(dwiNUM_OP+1)*(t+1)]/(double)(avg_TWI*Vol));                  
#endif                                                           
	  }
	}

	dwivec_clearer(avg_vector);      
	
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
