/*MASS TEST PROGRAM: it has been created in order to evaluate the renormalized mass of th theory */



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

  int i,k;
  int counter=0;
  int ncnfg;
  double DeltaT;
  double mean_mass=0.0;
  double Volume_mean=0.0;
  double err_mass=0.0;
  double accettanza, accettanza_cluster;
  char string[100];
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
  
  
  
  
  
  /* MEASUREMENT */
  
  
  
  mylog("[MAIN]: Starting with the measurements\n");
  
  double momp[3]={2*mPI/Nt,0.0,0.0};
  double mom0[3]={0.0,0.0,0.0};
  complex fp;
  complex f0;
  double fp2=0.0;
  double f02=0.0;
  double *massa=(double*)malloc(N_traj*sizeof(double));
  
  for(k=0;k<N_traj;k++){
    
    counter++;
    gettimeofday(&t1, NULL);

    for(i=0;i<N_meas;i++){
    
    update_sweep(N_sweeps,sweep_ratio);  
    phi_tilda(momp,Phi,&fp);   
    phi_tilda(mom0,Phi,&f0);
    fp2 += (fp.Re*fp.Re+fp.Im*fp.Im);
    f02 += (f0.Re*f0.Re);

    }

    fp2 /= N_meas;
    f02 /= N_meas;
    
    gettimeofday(&t2, NULL);
    
    DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
    massa[k]=sqrt(momp[0]*momp[0]*fp2/(f02-fp2));
    mean_mass+=massa[k]/N_traj;
    err_mass+=massa[k]*massa[k]/N_traj;
      
    mylog("<phi(p) square> %e <phi(0) square> %e \n",fp2,f02);
    mylog("Renormalized mass = %e generated in %e seconds\n",massa[k],(double)DeltaT/1000000.0);
    mylog("\n");
    
    fp2=0;
    f02=0;

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
    
    if(counter%N_save==0){
      
      Volume_mean=mean(Phi,Vol);
      ncnfg=(k+1)+N_term+Last_count;      
      
      mylog("[MAIN]: Saving field configuration\n");
      
      gettimeofday(&t1, NULL);
      sprintf(string,"%s/%dx%dx%dx%d_m%lf_l%lf_g%lf_k%lf_%d.bin","Cnfg",Nt,Nx,Ny,Nz,mass,lambda,g,kappa,ncnfg);
      write_output(string,Phi,Volume_mean);
      gettimeofday(&t2, NULL);
      
      DeltaT=(t2.tv_sec*1000000+t2.tv_usec)-(t1.tv_sec*1000000+t1.tv_usec);
      mylog("[MAIN]: Configuration %d saved in %lf seconds\n",ncnfg,(double)DeltaT/1000000.0);
      mylog("\n");
      
    }
  }
  
  gettimeofday(&T2, NULL);
  
  DeltaT=(T2.tv_sec*1000000+T2.tv_usec)-(T1.tv_sec*1000000+T1.tv_usec);
  err_mass=sqrt((err_mass-mean_mass*mean_mass)/(N_traj-1));
  
  mylog("[MAIN]: Renormalized mass = %lf +/- %lf \n",mean_mass,err_mass);
  mylog("[MAIN]: Total simulation time = %lf seconds \n",(double)DeltaT/1000000.0);
  
  
  
  
  
  /* FREEING MEMORY/CLOSING OUTPUT FILES*/ 
  
  
  
  mylog("[MAIN]: Freeing memory, closing files");
  mylog("\n");
  
  free_geo();
  free_cluster();
  free_Lattice(); 
  
  mylog("[MAIN]: closing program\n");
  
  return 0;
  
}
