/* INPUT_OUTPUT: library for input/output function definition and implementation*/



#include"input_output.h"
#include<stdlib.h>
#include<math.h>
#include"geometry.h"

/* INIT_MYLOG: function used to define the log file                               */
/* -take the logfile name as input                                                */
/* IF-the user does not give a name for the log file                              */
/*   -open STARTUP_ERROR file: a name need to be specified                        */
/* ELSE IF- the Cnfg directory does not exist                                     */
/*   -open STARTUP_ERROR file: Cnfg directory needs to be created                 */
/* ELSE IF- the Log directory does not exist                                      */
/*   -open STARTUP_ERROR file: Log directory needs to be created                  */
/* -open the logfile inside the directory Log                                     */
/* IF-the file does not exist: create it                                          */
/* ELSE-the file does exist                                                       */
/*   -open STARTUP_ERROR file: such a logfile already exist, cannot overwrite     */



void init_mylog(char *my_FILE,int narg){


  DIR* cnfg_dir=NULL;
  DIR* log_dir=NULL;
  FILE* streamtest;
  char string[100];

  if((cnfg_dir=opendir("Cnfg"))==NULL){
    
    stream = fopen ("STARTUP_ERROR","w");
    fprintf(stream,"[I/O] Cnfg directory does not exist, create it!\n");
    fprintf(stream,"[I/O] Run stopped \n");
    exit(1);
    
  }else if((log_dir=opendir("Log"))==NULL){
    
    stream = fopen ("STARTUP_ERROR","w");
    fprintf(stream,"[I/O] Log directory does not exist, create it!\n");
    fprintf(stream,"[I/O] Run stopped \n");
    exit(1);
    
  }
  
  
  if(narg<2){
    
    sprintf(string,"%s","Log/std_logfile");
    streamtest = fopen (string,"r");
    
    if(!streamtest){
      
      stream = fopen(string,"w");
      fprintf(stream,"[I/O] STANDARD LOG FILE\n");  
      
    }else{
      
      fclose(streamtest);
      closedir(log_dir);
      stream = fopen ("STARTUP_ERROR","w");
      fprintf(stream,"[I/O] Standard log file already exist!\n");
      fprintf(stream,"[I/O] Cancel std_logfile o give a new file name\n");
      fprintf(stream,"[I/O] Run stopped \n");
      fclose(stream);
      exit(1);
      
    }    
  }else{ 
    
    sprintf(string,"%s/%s","Log",my_FILE);  
    streamtest=fopen(string,"r");
    
    if(!streamtest){
      
    stream = fopen(string,"w");
    fprintf(stream,"[I/O] NEW LOG FILE\n");  
    
    }else{
      
      fclose(streamtest);
      closedir(log_dir);
      stream = fopen ("STARTUP_ERROR","w");
      fprintf(stream,"[I/O] Log file already exist!\n");
      fprintf(stream,"[I/O] Run stopped \n");
      fclose(stream);
      exit(1);
      
    }    
  }
}




/* MYLOG: function used to print report informations inside the logfile                              */
/* -define a customized fprintf that send the report infos into the logfile, using variadic function */



void mylog( const char * report_info, ... ){

  va_list report_parmters;
  va_start(report_parmters,report_info);

  vfprintf(stream,report_info,report_parmters);

  va_end(report_parmters);
  fflush(stream);

}





/* SEARCHWORD */



int searchword(FILE* SourceFile,char* searchingString){
   
  char LineInput[512];
  int sizeLine,stringsize;
  
  rewind(SourceFile);
  stringsize = strlen(searchingString);
  while (fgets(LineInput, 512, SourceFile)!=NULL)
    {
      
      if(strncmp (LineInput,searchingString,stringsize)==0){
	sizeLine = strlen(LineInput);
	fseek(SourceFile,-sizeLine+stringsize,SEEK_CUR);
	return true;
	
      }
    }
  
  mylog("[I/O]: Voice %s not found inside parameters file\n",searchingString);
  mylog("[I/O]: Run stopped \n");
  mylog("\n");
  return false;
  exit(1);
  
}





/* ERROR: function used to printf error statement during fiel reading */
/* IF-the value of bool is true (e.g. a parameter value is missing)   */
/*   -print information about missing value (given in statString)     */
/* -stop the run                                                      */
/* -report to logfile                                                 */



void error(bool stat, char* statString){
  
  if(stat==true){
    
    mylog(statString);
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }
}





/* READINPUT                                            */
/* -take filename, check it existence                   */
/* -open file                                           */
/* -read the following input values                     */
/*  1-dimension of the system                           */
/*  2-lattice extension                                 */
/*  3-field dimension                                   */
/*  4-number of measurement                             */
/*  5-number of termalization steps                     */
/*  6-number of sweeps steps                            */
/*  7-mass parameter                                    */
/*  8-lambda parameter                                  */
/*  9-g parameter                                       */
/*  10-kappa parameter                                  */
/*  11-metropolis control parameter                     */
/*  12-integration interval of Wilson Flow              */
/*  13-number of Wilson Flow integration steps          */
/*  14-rate measurement during flow evolution           */
/*  15-status variable(used for setting initial fields) */
/* -evaluate system volume for different dimensions     */
/* -close file                                          */



void readinput(char *my_FILE){

  FILE *fp_my_file;
  status=(char*)malloc(100*sizeof(char)); 
  mass2_sign=(char*)malloc(10*sizeof(char)); 

  fp_my_file = fopen(my_FILE,"r");

  if(!fp_my_file){

    mylog("[I/O]: Cannot find parameters file!\n");
    mylog("[I/O]: Run stopped\n");
    exit(1);

  }else{

    mylog("[I/O]: Parameters filename ok\n");
    mylog("\n");
  }

  searchword(fp_my_file,"system dimensions = ");
  error(fscanf(fp_my_file,"%d \n",&d)!=1,"[I/O]: Missing system dimensions value\n"); 
 
  searchword(fp_my_file,"Nt = ");
  error(fscanf(fp_my_file,"%d \n",&Nt)!=1,"[I/O]: Missing Nt value\n");

  searchword(fp_my_file,"Nx = ");
  error(fscanf(fp_my_file,"%d \n",&Nx)!=1,"[I/O]: Missing Nx value\n");

  searchword(fp_my_file,"Ny = ");
  error(fscanf(fp_my_file,"%d \n",&Ny)!=1,"[I/O]: Missing Ny value\n");

  searchword(fp_my_file,"Nz = ");
  error(fscanf(fp_my_file,"%d \n",&Nz)!=1,"[I/O]: Missing Nz value\n");

  searchword(fp_my_file,"field dimension = ");
  error(fscanf(fp_my_file,"%d \n",&d_field)!=1,"[I/O]: Missing field dimension value \n");                 
  
  searchword(fp_my_file,"trajectories = ");
  error(fscanf(fp_my_file,"%d \n",&N_traj)!=1,"[I/O]: Missing number of trajectories \n");                   
  
  searchword(fp_my_file,"thermalization steps = ");
  error(fscanf(fp_my_file,"%d \n",&N_term)!=1,"[I/O]: Missing number of thermalization steps\n");                                      
  
  searchword(fp_my_file,"sweeps = ");
  error(fscanf(fp_my_file,"%d \n",&N_sweeps)!=1,"[I/O]: Missing number of sweeps steps\n");                   

  searchword(fp_my_file,"sweep_ratio = ");
  error(fscanf(fp_my_file,"%d \n",&sweep_ratio)!=1,"[I/O]: Missing sweep ratio value \n");                   

  searchword(fp_my_file,"report rate = ");
  error(fscanf(fp_my_file,"%d \n",&N_report)!=1,"[I/O]: Missing report rate value \n");                   

  searchword(fp_my_file,"saving rate = ");
  error(fscanf(fp_my_file,"%d \n",&N_save)!=1,"[I/O]: Missing saving rate value \n");                   
  
  searchword(fp_my_file,"measure rate = ");
  error(fscanf(fp_my_file,"%d \n",&N_meas)!=1,"[I/O]: Missing measure rate value \n");                   
  
  searchword(fp_my_file,"mass = ");
  error(fscanf(fp_my_file,"%lf \n",&mass)!=1,"[I/O]: Missing mass value \n");              
  
  searchword(fp_my_file,"lambda = ");
  error(fscanf(fp_my_file,"%lf \n",&lambda)!=1,"[I/O]: Missing lambda coupling value \n");                 
  
  searchword(fp_my_file,"g = ");
  error(fscanf(fp_my_file,"%lf \n",&g)!=1,"[I/O]: Missing g coupling value \n");                 
  
  searchword(fp_my_file,"kappa = ");
  error(fscanf(fp_my_file,"%lf \n",&kappa)!=1,"[I/O]: Missing kappa coupling value \n");                 
 
  searchword(fp_my_file,"square mass sign = ");  
  error(fscanf(fp_my_file,"%s \n",mass2_sign)==0,"[I/O]: Missing square mass sign");
  
  searchword(fp_my_file,"Integrator order = ");
  error(fscanf(fp_my_file,"%d \n",&int_ord)!=1,"[I/O]: Missing integrator order \n");
  
  searchword(fp_my_file,"Bulk mass = ");
  error(fscanf(fp_my_file,"%lf \n",&Bulk_mass)!=1,"[I/O]: Missing Bulk mass value \n");                   
  
  searchword(fp_my_file,"Bulk lambda = ");
  error(fscanf(fp_my_file,"%lf \n",&Bulk_lambda)!=1,"[I/O]: Missing Bulk lambda coupling value \n");                 
  
  searchword(fp_my_file,"Bulk g = ");
  error(fscanf(fp_my_file,"%lf \n",&Bulk_g)!=1,"[I/O]: Missing Bulk g coupling value \n");                 
  
  searchword(fp_my_file,"Bulk kappa = ");
  error(fscanf(fp_my_file,"%lf \n",&Bulk_kappa)!=1,"[I/O]: Missing Bulk kappa coupling value \n");                 
  
  searchword(fp_my_file,"delta = ");
  error(fscanf(fp_my_file,"%lf \n",&delta)!=1,"[I/O]: Missing Metropolis delta value \n");                  

  searchword(fp_my_file,"random seed = ");
  error(fscanf(fp_my_file,"%d \n",&seed)!=1,"[I/O]: Missing random seed value \n");                  
  
  searchword(fp_my_file,"Wilson Flow dt = ");
  error(fscanf(fp_my_file,"%lf \n",&dt)!=1,"[I/O]: Missing Wilson Flow integration dt \n");                  
  
  searchword(fp_my_file,"Wilson Flow steps = ");
  error(fscanf(fp_my_file,"%d \n",&N_steps)!=1,"[I/O]: Missing Wilson Flow integration steps \n");                  

  searchword(fp_my_file,"Wilson Flow rate = ");
  error(fscanf(fp_my_file,"%d \n",&N_flow)!=1,"[I/O]: Missing Wilson Flow measure rate \n");                  
  
  searchword(fp_my_file,"start = ");
  error(fscanf(fp_my_file,"%s",status)==0,"[I/O]: Missing start condition \n");
  
  fclose(fp_my_file); 
  

  if(d==2){
  
    Nz=0;
    Ny=0;
    Vol=Nt*Nx;

  }else if(d==3){
  
    Nz=0;
    Vol=Nt*Nx*Ny;
    
  }else if(d==4){
    
    Vol=Nt*Nx*Ny*Nz;
    
  }

  
  if(strcmp(mass2_sign,"+")==0){
    
    mass_square=mass*mass;
    
  }else if(strcmp(mass2_sign,"-")==0){
    
    mass_square=-mass*mass;
    
  }else if(strcmp(mass2_sign,"+")!=0||strcmp(mass2_sign,"-")!=0){
    
    mylog("[I/O]: You need to specify a sign for the square mass \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }  

  if(int_ord<1){

    mylog("[I/O]: You need to specify which kind of integrator you want to use for the Wilson Flow \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);

  }
  
  if(N_report==0||N_save==0||N_meas==0||N_report>N_term){
    
    mylog("[I/O]: Report rate, save rate, measure rate can't be zero or greater than the total number of trajectories \n");
    mylog("[I/O]: Number of trajectories, N_traj = %d\n",N_traj);
    mylog("[I/O]: Report rate, N_report = %d\n",N_report);
    mylog("[I/O]: Report rate, N_save = %d\n",N_save);
    mylog("[I/O]: Report rate, N_meas = %d\n",N_meas);
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }  
}


/* WRITE_OUTPUT */



void write_output(char *my_FILE, double *Obs, double Obs_check){

  int out1,out2,out3,out4,out5,out6;
  FILE* fp_my_File;

  fp_my_File=fopen(my_FILE,"wb");
  out1=fwrite(&Obs_check,sizeof(double),1,fp_my_File);
  out2=fwrite(&mass,sizeof(double),1,fp_my_File);
  out3=fwrite(&lambda,sizeof(double),1,fp_my_File);
  out4=fwrite(&g,sizeof(double),1,fp_my_File);
  out5=fwrite(&kappa,sizeof(double),1,fp_my_File);
  out6=fwrite(Obs,sizeof(double),Vol,fp_my_File);
  
  if(out6!=Vol||out1!=1||out2!=1||out3!=1||out4!=1||out5!=1){
    
    mylog("[I/O]: Error during output writing\n");
    mylog("\n");
    exit(1);
    
  }
  
  fclose(fp_my_File);
  
}





/* FLOWINPUT */
/* -take filename, check it existence                   */
/* -open file                                           */
/* -read the following input values                     */
/* 1-tipe of flow time unity measure (t or c(t))        */
/* 2-number of flow times                               */
/* 3flow times                                          */



#define rounding(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))




void flow_input(char*my_FILE){

  FILE *fp_my_file;    
  int count;
  double value;
  
  char* flowdef=(char*)malloc(10*sizeof(char)); 
  
  
  fp_my_file = fopen(my_FILE,"r");
  
  if(!fp_my_file){
    
    mylog("[I/O]: Cannot find flow times file!\n");
    mylog("[I/O]: Run stopped\n");
    exit(1);

  }else{

    mylog("[I/O]: Flowtimes filename ok\n");
    mylog("\n");
  }

  searchword(fp_my_file,"flow time definition = ");
  error(fscanf(fp_my_file,"%s",flowdef)==0,"[I/O]: Missing relative flow times definition \n");

  searchword(fp_my_file,"number of flow times = ");
  error(fscanf(fp_my_file,"%d",&N_flowtimes)!=1,"[I/O]: Missing number of flow times \n");  
  Flow_vec=(int*)malloc(N_flowtimes*sizeof(int));

  
  if(strcmp(flowdef,"t")==0){

    count=0;
 
    while(fscanf(fp_my_file,"%lf\n",&value)==1){
      
      Flow_vec[count]=rounding(value/dt);
      count++;
      
    }    
  }else if(strcmp(flowdef,"c(t)")==0){

    count=0;    
    
    while(fscanf(fp_my_file,"%lf\n",&value)==1){
      
      Flow_vec[count]=rounding((value*value*(double)(Nt*Nt)/8.)/dt);
      count++;
      
    }        
  }else{
    
    mylog("[I/O]: Incorrect relative flow times definition \n");
    mylog("[I/O]: Available options are  \n");
    mylog("[I/O]: (1) 't' take flow times as given by the input file  \n");
    mylog("[I/O]: (2) 'c(t)' take  (smearing radius)/(lattice size) ratios and convert them in flow times \n");
    mylog("[I/O]: Run stopped\n");
    exit(1);
    
  }
  
  if (count!=N_flowtimes){
    
    mylog("[I/O]: The number of given flow times is not correct \n");
    mylog("[I/O]: counter = %d, N_flowtimes = %d \n",count, N_flowtimes);
    mylog("[I/O]: Run stopped\n");
    exit(1);
    
  }
  
  if(Flow_vec[0]==0){
    
    mylog("[I/O]: Your dt is too large, cannot assign a non zero value to %lf \n", Flow_vec[0]*dt);
    mylog("[I/O]: Change dt or the series of flow times \n");
    mylog("[I/O]: Run stopped\n");
    exit(1);
    
  }
  
  int check=Flow_vec[0];
  
  for(count=1;count<N_flowtimes;count++){
    
    if(Flow_vec[count]==check){
      
      mylog("[I/O]: Your dt is too large, cannot distinguish between two flowtimes\n");
      mylog("[I/O]: Flow_vec[%d] = %lf ; Flow_vec[%d] = %lf \n",count-1,check*dt,count,check*dt);      
      mylog("[I/O]: Change dt or the series of flow times \n");
      mylog("[I/O]: Run stopped\n");
      exit(1);
      
    }

    check=Flow_vec[count];
    
  }
  
  N_steps=Flow_vec[N_flowtimes-1];

  mylog("[I/O]: Number of flow steps modified to cover all the flow times from input \n");
  mylog("[I/O]: Nsteps = %d \n",N_steps);


  fclose(fp_my_file); 
  
}






/* TWI_INPUT */
/* -take filename, check it existence                   */
/* -open file                                           */
/* -read the following input values                     */
/* 1-point where we apply the translation               */
/* 2-probe observables position                         */
/* LHS measure options                                  */
/* RHS measure options                                  */


void TDWI_input(char *my_FILE){

  FILE *fp_my_file;    
  char *lhs_m=(char*)malloc(10*sizeof(char));
  char *rhs_m=(char*)malloc(10*sizeof(char));
  
  fp_my_file = fopen(my_FILE,"r");
  if(!fp_my_file){
   
    mylog("[I/O]: Cannot find TWI file!\n");
    mylog("[I/O]: Run stopped\n");
    exit(1);
    
  }else{
    
    mylog("[I/O]: TWI filename ok\n");
    mylog("\n");
  }
  
  int dist=0;
  searchword(fp_my_file,"number of distances = ");
  error(fscanf(fp_my_file,"%d",&dist)==0,"[I/O]: Missing number of distances \n");
  if(dist<=0){
    mylog("[I/O]: Number of distances must be positive \n");
    mylog("[I/O]: Run stopped\n");
    exit(1);
  }
  
  ndistances=dist;
  
  searchword(fp_my_file,"rhs measure = ");
  error(fscanf(fp_my_file,"%s",rhs_m)==0,"[I/O]: Missing rhs measuring condition (on/off) \n");
  
  searchword(fp_my_file,"lhs measure = ");
  error(fscanf(fp_my_file,"%s",lhs_m)==0,"[I/O]: Missing lhs measuring condition (on/off) \n");

  searchword(fp_my_file,"average block = ");
  error(fscanf(fp_my_file,"%d",&avg_TWI)==0,"[I/O]: Missing block size for intermediate averages \n");
  
  if(strcmp(lhs_m,"on")==0){
    
    lhs_measure=1;
    
  }else if(strcmp(lhs_m,"off")==0){
    
    lhs_measure=0;
    
  }else{
    
    mylog("[I/O]: You need to specify if you want to measure the LHS (on/off) \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }  
  
  if(strcmp(rhs_m,"on")==0){
    
    rhs_measure=1;
    
  }else if(strcmp(rhs_m,"off")==0){
    
    rhs_measure=0;
    
  }else{
    
    mylog("[I/O]: You need to specify if you want to measure the RHS (on/off) \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }

  if(strcmp(rhs_m,"off")==0 && strcmp(lhs_m,"off")==0 ){
    
    mylog("[I/O]: You need to measure something idiot \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }
  
  
  
  if(avg_TWI<=0||N_meas%avg_TWI!=0){
    
    mylog("[I/O]: Block size has to be a positive divisor of the number of measurements \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }  
  
}




void TWI_input(char*my_FILE){

  FILE *fp_my_file;    
  int traslt,traslx,trasly,probet,probex,probey;
  char *lhs_m=(char*)malloc(10*sizeof(char));
  char *rhs_m=(char*)malloc(10*sizeof(char));
  
  fp_my_file = fopen(my_FILE,"r");
  
  if(!fp_my_file){
    
    mylog("[I/O]: Cannot find TWI file!\n");
    mylog("[I/O]: Run stopped\n");
    exit(1);

  }else{

    mylog("[I/O]: TWI filename ok\n");
    mylog("\n");
  }


  

  searchword(fp_my_file,"translation coord t = ");
  error(fscanf(fp_my_file,"%d",&traslt)==0,"[I/O]: Missing translation t coordinate \n");

  searchword(fp_my_file,"translation coord x = ");
  error(fscanf(fp_my_file,"%d",&traslx)==0,"[I/O]: Missing translation x coordinate \n");

  searchword(fp_my_file,"translation coord y = ");
  error(fscanf(fp_my_file,"%d",&trasly)==0,"[I/O]: Missing translation y coordinate \n");

  searchword(fp_my_file,"probe coord t = ");
  error(fscanf(fp_my_file,"%d",&probet)==0,"[I/O]: Missing probe t coordinate \n");

  searchword(fp_my_file,"probe coord x = ");
  error(fscanf(fp_my_file,"%d",&probex)==0,"[I/O]: Missing probe x coordinate \n");

  searchword(fp_my_file,"probe coord y = ");
  error(fscanf(fp_my_file,"%d",&probey)==0,"[I/O]: Missing probe y coordinate \n");
  
  if(traslt<0||traslx<0||trasly<0||probet<0||probex<0||probey<0){
    
    mylog("[I/O]: Probe and translation coordinates must be positive \n");
    mylog("[I/O]: Run stopped\n");
    exit(1);
  }



  if(traslt>=Nt||traslx>=Nx||trasly>=Ny||probet>=Nt||probex>=Nx||probey>=Ny){

    mylog("[I/O]: Probe and translation coordinates must be less than lattice extension \n");
    mylog("[I/O]: Run stopped\n");
    exit(1);
  }



  traslpos=get_index(traslt,traslx,trasly);
  probepos=get_index(probet,probex,probey);          


  searchword(fp_my_file,"rhs measure = ");
  error(fscanf(fp_my_file,"%s",rhs_m)==0,"[I/O]: Missing rhs measuring condition (on/off) \n");

  searchword(fp_my_file,"lhs measure = ");
  error(fscanf(fp_my_file,"%s",lhs_m)==0,"[I/O]: Missing lhs measuring condition (on/off) \n");

  searchword(fp_my_file,"average block = ");
  error(fscanf(fp_my_file,"%d",&avg_TWI)==0,"[I/O]: Missing block size for intermediate averages \n");





  if(strcmp(lhs_m,"on")==0){
    
    lhs_measure=1;
    
  }else if(strcmp(lhs_m,"off")==0){
    
    lhs_measure=0;
    
  }else{
    
    mylog("[I/O]: You need to specify if you want to measure the LHS (on/off) \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }  

  if(strcmp(rhs_m,"on")==0){
    
    rhs_measure=1;
    
  }else if(strcmp(rhs_m,"off")==0){
    
    rhs_measure=0;
    
  }else{
    
    mylog("[I/O]: You need to specify if you want to measure the RHS (on/off) \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);
    
  }

  if(strcmp(rhs_m,"off")==0 && strcmp(lhs_m,"off")==0 ){

    mylog("[I/O]: You need to measure something idiot \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);

  }



  if(avg_TWI<=0||N_meas%avg_TWI!=0){

    mylog("[I/O]: Block size has to be a positive divisor of the number of measurements \n");
    mylog("[I/O]: Run stopped \n");
    mylog("\n");
    exit(1);

}  
  
  return;
}



/* SEED_READINPUT                                                         */
/* -read the value of the current random seed from the file "random_seed" */



void seed_readinput(char *my_FILE){

  int out;
  FILE *fp_my_file; 
  fp_my_file = fopen(my_FILE,"r");

  if(!fp_my_file){

    perror("Idiot\n");

  }else{

    printf("random seed filename ok\n");

  }

 out=fscanf(fp_my_file,"%d",&seed);
 if(out!=1) exit(1);      
 fclose(fp_my_file);

}





/* SEED_OUTPUT                                             */
/* -save the updated random seed in the file "random_seed" */



void seed_output(char *my_FILE){
  
  int out;
  FILE *fp_my_file; 
  fp_my_file = fopen(my_FILE,"w");
  
  if(!fp_my_file){
    
    perror("Idiot\n");
    
  }else{
    
    printf("random seed filename ok\n");
    
  }
  
  out=fprintf(fp_my_file,"%d",seed+1);   
  if(out==0) exit(1);                 
  fclose(fp_my_file);
  
}





