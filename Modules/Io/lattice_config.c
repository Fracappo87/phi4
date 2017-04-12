/*LATTICE_CONFIG: configures the initial values of fields in every lattice site*/



#include"lattice_config.h"
#include"ranlux.h"
#include"input_output.h"
#include<stdio.h>
#include<stdlib.h>
//#include<malloc.h>
#include<string.h>



 
/* CONFIGURATION                                                                                 */
/* -allocate memory for field values vector                                                      */
/* IF-variable status="start=cold", then set cold start for the simulation                       */
/*   L-for every lattice site                                                                    */ 
/*    -set the values of the field at 0.5                                                        */
/* IF-variable status="start=hot", then set hot start for the simulation                         */
/*   -generate a one component random vector rand                                                */
/*   L-for every lattice site                                                                    */
/*    -set the values of the field equal to rand[0]                                              */
/* IF-variable status="start=mass_lambda_..._volume.bin", then set hot start for the simulation  */
/*   -open file with name status                                                                 */
/*   L-for every lattice site                                                                    */
/*    -initialize the field values                                                               */
/*   -close file                                                                                 */



void configuration(){

  int i;
  double check_value1,check_value2,check_value3,check_value4,check_value5;
  double Volume_mean=0;
  char string[100];
  Phi=(double*)malloc(Vol*sizeof(double));

  mylog("[LATTICE_CONFIG]: Setting lattice field configuration\n");
  mylog("[LATTICE_CONFIG]: Number of field components, d = %d\n",d_field);
  
  if(strcmp(status,"cold")==0){ 
    
    mylog("[LATTICE_CONFIG]: Cold start\n");
    mylog("\n");
 
    for(i=0;i<Vol;i++){    
      
   iPhi(i)=0.5;
   
    }
  }else if(strcmp(status,"hot")==0){
    
    double rand[1];
     mylog("[LATTICE_CONFIG]: Hot start\n");
     mylog("\n");

    for(i=0;i<Vol;i++){    
      
      ranlxd(rand,1);
      iPhi(i)=2.*rand[0]-1.;
       
      
    }
  }else if((strcmp(status,"hot")!=0)||(strcmp(status,"cold")!=0)){
    
    FILE* fp;
    int out1,out2,out3,out4,out5,out6;
    int check1,check2,check3,check4;
    double check5,check6,check7,check8;

    sprintf(string,"%s/%s","Cnfg",status);
    sscanf(status,"%dx%dx%dx%d_m%lf_l%lf_g%lf_k%lf_%d.bin",&check1,&check2,&check3,&check4,&check5,&check6,&check7,&check8,&Last_count);
    fp=fopen(string,"rb");

    if(check1!=Nt||check2!=Nx||check3!=Ny||check4!=Nz||check5!=mass||check6!=lambda||check7!=g||check8!=kappa){

	mylog("[LATTICE_CONFIG]: You are not starting from a configuration with the same parameters, idiot\n");
	mylog("[LATTICE_CONFIG]: Nt  %d %d\n",check1, Nt);
	mylog("[LATTICE_CONFIG]: Nx  %d %d\n",check2, Nx);
	mylog("[LATTICE_CONFIG]: Ny  %d %d\n",check3, Ny);
	mylog("[LATTICE_CONFIG]: Nz  %d %d\n",check4, Nz);
	mylog("[LATTICE_CONFIG]: Mass parameter %lf %lf\n",check5, mass);
	mylog("[LATTICE_CONFIG]: Lambda couplig %lf %lf\n",check6, lambda);
	mylog("[LATTICE_CONFIG]: G coupling %lf %lf\n",check7, g);
	mylog("[LATTICE_CONFIG]: Kappa coupling %lf %lf\n",check8, kappa);	
	mylog("[LATTICE_CONFIG]: Run stopped ");
	exit(1);

    }
    
    mylog("[LATTICE_CONFIG]: Start from saved configuration\n");
        
    if(!fp){
      
      mylog("[LATTICE_CONFIG]: Cannot find configuration file!\n");
      mylog("[LATTICE_CONFIG]: Possible start options are: \n");
      mylog("[LATTICE_CONFIG]: start=hot \n");
      mylog("[LATTICE_CONFIG]: start=cold \n");
      mylog("[LATTICE_CONFIG]: start=filename \n");
      mylog("[LATTICE_CONFIG]: Run stopped \n");
      exit(1);
      
    }else{
      
      mylog("[LATTICE_CONFIG]: Initial configuration taken from %s\n",string);
      
      out1=fread(&check_value1,sizeof(double),1,fp);
      out2=fread(&check_value2,sizeof(double),1,fp);
      out3=fread(&check_value3,sizeof(double),1,fp);
      out4=fread(&check_value4,sizeof(double),1,fp);
      out5=fread(&check_value5,sizeof(double),1,fp);
      out6=fread(Phi,sizeof(double),Vol,fp);

      for(i=0;i<Vol;i++) Volume_mean+=iPhi(i)/Vol;

      mylog("[LATTICE_CONFIG]: Checking parameters\n");
      mylog("\n");

      if(out6!=Vol||out1!=1||out2!=1||out3!=1||out4!=1||out5!=1||check_value1!=Volume_mean||check_value2!=mass||check_value3!=lambda||check_value4!=g||check_value5!=kappa){
	
	mylog("[LATTICE_CONFIG]: You are not starting from a configuration with the same parameters, idiot\n");
	mylog("[LATTICE_CONFIG]: Field volume mean %lf %lf\n",check_value1, Volume_mean);
	mylog("[LATTICE_CONFIG]: Mass parameter %lf %lf\n",check_value2, mass);
	mylog("[LATTICE_CONFIG]: Lambda couplig %lf %lf\n",check_value3, lambda);
	mylog("[LATTICE_CONFIG]: G coupling %lf %lf\n",check_value4, g);
	mylog("[LATTICE_CONFIG]: Kappa coupling %lf %lf\n",check_value5, kappa);	
	mylog("[LATTICE_CONFIG]: Run stopped ");
	exit(1);
      
      }

      fclose(fp);
      
    }
  }
}



/*FREE_LATTICE                                    */
/* -deallocate the memory for field values vector */



void free_Lattice(){
  
  free(Phi);
  Phi=NULL;

  mylog("[LATTICE]: field vector's memory freed\n");
  mylog("\n");

}
