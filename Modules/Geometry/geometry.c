/*GEOMETRY: library for system geometry initialization */
/*-coordinate convention (0-t),(1-x),(2-y),(3-z)       */



#include<stdio.h>
#include<math.h>
#include"geometry.h"




/*COORDINATES                                             */
/* -extract the coordinates values from the lattice index */





void coordinates(int i,int* x){
  
  x[0]= i%Nt;
  x[1]= (i/Nt)%Nx;
  x[2]= (i/(Nt*Nx))%Ny;
  
}





/*GET_INDEX                                         */
/* - assign an integer index for every lattice site */



int get_index(int t,  int x, int y){

  return t+Nt*x+Nt*Nx*y;
   
}




/* GEOMETRY ALGORITHM                                                                   */
/* -allocate memory for next_nbr coordinate vectors d_up, d_down                        */
/* -macro: substitute iup(site,dir) and idn(site,dir) with the following slot element   */
/* L-for every t                                                                        */
/*  L-for every x                                                                       */
/*   L-for every y                                                                      */
/*    -define translation unit vector dm[mu]                                            */
/*    -evaluate the coordnates of next neighbour points with respect the point (t,x,y)  */
/*    -impose P.B.C using %                                                             */
/*    -get the index of the next-nbr point, allocate it as the (i*d+/-mu)th component   */
/* L-for every y                                                                        */
/*  -tabulate the structure coordvec with the entries of the vectors that label all     */
/*   the lattice points                                                                 */
/* -report to logfile                                                                   */



void geometry(){

  int t,x,y;
  int mu;
  d_up=(int*)malloc(d*Vol*sizeof(int));
  d_down=(int*)malloc(d*Vol*sizeof(int));

  coordvec=(coord*)malloc(Vol*sizeof(coord));

  
#define iup(site,dir) d_up[(site)*d+(dir)] 
#define idn(site,dir) d_down[(site)*d+(dir)]

  mylog("[GEOMETRY]: Setting system geometry with PBC\n");
 
  for(t=0;t<Nt;t++)
    for(x=0;x<Nx;x++)
      for(y=0;y<Ny;y++)
	for(mu=0;mu<d;mu++){
	
	  int dmu[3]={0,0,0};
	 
	  dmu[mu]=1;
	  int tn = (t+dmu[0])%Nt;
	  int xn = (x+dmu[1])%Nx;
	  int yn = (y+dmu[2])%Ny;
	  iup(get_index(t,x,y),mu)=get_index(tn,xn,yn);

	  dmu[mu]=-1;
	  tn = (t+dmu[0]+Nt)%Nt;
	  xn = (x+dmu[1]+Nx)%Nx;
	  yn = (y+dmu[2]+Ny)%Ny;
	  idn(get_index(t,x,y),mu)=get_index(tn,xn,yn);
  
	}  

  mylog("[GEOMETRY]: System dimension = %d\n",d);
  mylog("[GEOMETRY]: Nt = %d, Nx = %d, Ny = %d, Nz = %d\n",Nt,Nx,Ny,Nz);
  mylog("[GEOMETRY]: Volume = %d\n",Vol);

  mylog("[GEOMETRY]: Filling vector coordvec");
  mylog("\n");
  for(t=0;t<Vol;t++)
    coordinates(t,coordvec[t].x);

}




/*FREE_GEO                                   */
/* -deallocate the next-nb vectors' memory   */



void free_geo(){

  free(d_down);
  free(d_up);
  d_up=NULL;
  d_down=NULL;
  
  mylog("[GEOMETRY]: nnb vector's memory freed\n");
  mylog("\n");

}
