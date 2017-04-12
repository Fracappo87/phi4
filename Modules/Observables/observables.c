/* OBSERVABLES: library for the definitions of measured observables */



#include <math.h>
#include <stdio.h>
#include "observables.h"
#include "geometry.h"
#include "lattice_config.h"





/* MEAN: define the mean of a lattice quantity */



double mean(double* vec, int d){  
  int i;
  double sum=0;
  
  for(i=0;i<d;i++) sum+=vec[i];

  return sum/(double)d;

}





/* CORR: function for the evaluation of correlators                                                                        */
/* In this case, since we can use the translational invariance of the given geometry,                                      */
/* we can evaluate the volume mean of the correlation function, once the distance "dis" between the operators is given     */
/* L-for an index i lesser than the correlation distance                                                                   */
/*  -define the correlation distance along the different dimension, iteratively constructed using the nnb vectors          */
/*  L-for each point along Nt                                                                                              */
/*   L-for each point along Nx                                                                                             */
/*    L-for each point along Ny                                                                                            */
/*     -sum the correlators with respect the distance "dis" along the 3 different system dimensions                        */
/*    -make a rigid shift along the y axis of all the application points                                                   */
/*   -make a rigid shift along the x axis of all the application points                                                    */
/* Divide the result of the above sum by a volume factor                                                                   */



double corr(int dis, double *phi){
  
  int i,j,k;
  int pos=0,pos1=0,pos2=0,pos3=0;
  for(i=0;i<dis;i++){
    pos1=d_up[pos1*d];
    pos2=d_up[pos2*d+1];
    pos3=d_up[pos3*d+2];
  }
  double corr=0.;
  for(i=0;i<Nt;i++){
    for(j=0;j<Nx;j++){
      for(k=0;k<Ny;k++){
	
	corr+=Obs(pos)*(Obs(pos1)+Obs(pos2)+Obs(pos3));
	pos=d_up[pos*d+2];
	pos1=d_up[pos1*d+2];
	pos2=d_up[pos2*d+2];
	pos3=d_up[pos3*d+2];
      }
      pos=d_up[pos*d+1];
      pos1=d_up[pos1*d+1];
      pos2=d_up[pos2*d+1];
      pos3=d_up[pos3*d+1];
    }
    pos=d_up[pos*d];
    pos1=d_up[pos1*d];
    pos2=d_up[pos2*d];
    pos3=d_up[pos3*d];
  }
  return corr/(double)(Vol*3.);
  
}


double rhschainnew(int dis, double *phi){
  
  int i,j,k;
  int pos=0,pos1=0,pos2=0,pos3=0,pos4,pos5,pos6;
  for(i=0;i<dis;i++){
    pos1=d_up[pos1*d];
    pos2=d_up[pos2*d+1];
    pos3=d_up[pos3*d+2];
  }

  pos4=d_up[pos1*d];
  pos5=d_up[pos2*d+1];
  pos6=d_up[pos3*d+2];

  double corr=0.;
  for(i=0;i<Nt;i++){
    for(j=0;j<Nx;j++){
      for(k=0;k<Ny;k++){
	
	corr+=Obs(pos)*(DObs(pos1)*(phi[pos1]-phi[pos4])+DObs(pos2)*(phi[pos2]-phi[pos5])+DObs(pos3)*(phi[pos3]-phi[pos6]));
	pos=d_up[pos*d+2];
	pos1=d_up[pos1*d+2];
	pos2=d_up[pos2*d+2];
	pos3=d_up[pos3*d+2];
	pos4=d_up[pos4*d+2];
	pos5=d_up[pos5*d+2];
	pos6=d_up[pos6*d+2];
      }
      pos=d_up[pos*d+1];
      pos1=d_up[pos1*d+1];
      pos2=d_up[pos2*d+1];
      pos3=d_up[pos3*d+1];
      pos4=d_up[pos4*d+1];
      pos5=d_up[pos5*d+1];
      pos6=d_up[pos6*d+1];
    }
    pos=d_up[pos*d];
    pos1=d_up[pos1*d];
    pos2=d_up[pos2*d];
    pos3=d_up[pos3*d];
    pos4=d_up[pos4*d];
    pos5=d_up[pos5*d];
    pos6=d_up[pos6*d];
  }
  return corr/(double)(Vol*3.);
  
}


/* SPATIAL_MEAN: d-1 dimensional mean, used for TWI */
/* This function will give you the d-1 dimensional mean of a given observable (defined in the macro in observables.h), for a given point along the time direction */
/* Define two lattice indices pos1, pos2 related to the temporal coordinate of the observable                                                                     */
/* L-for every point along the x axis                                                                                                                             */
/*  -set pos1=pos2                                                                                                                                                */
/*  L-for every point along the y axis                                                                                                                            */
/*   -sum the observable in the given point pos1                                                                                                                  */
/*   -make shift the observable in the nnb of pos1 along the diretion y                                                                                           */
/*  -make shift the observable in the nnb of pos1 along the diretion x                                                                                            */
/* -Divide the result of the sum above by an area factor                                                                                                          */



double spatial_mean(int coord_0, double *phi){
  
  int i,j;
  int pos1=get_index(coord_0,0,0);
  int pos2=get_index(coord_0,0,0);
  double obsr=0.;
  
  for(i=0;i<Nx;i++){
    
    pos1=pos2;
    for(j=0;j<Ny;j++){
      
      obsr+=Obs(pos1);
      pos1=d_up[pos1*d+2];
    }
    pos2=d_up[pos2*d+1];
  }
  
  return obsr/((double)(Nx*Ny));
  
}









/* SPATIAL_MEAN_MUL: d-1 dimensional mean, used for TWI                                                                                                           */
/* This function will give you the d-1 dimensional mean of a product of observables, for a given point along the time direction                                   */
/* Define two lattice indices pos1, pos2 related to the temporal coordinate of the observable                                                                     */
/* L-for every point along the x axis                                                                                                                             */
/*  -set pos1=pos2                                                                                                                                                */
/*  L-for every point along the y axis                                                                                                                            */
/*   -sum the observable in the given point pos1                                                                                                                  */
/*   -make shift the observable in the nnb of pos1 along the diretion y                                                                                           */
/*  -make shift the observable in the nnb of pos1 along the diretion x                                                                                            */
/* -Divide the result of the sum above by an area factor                                                                                                          */



double spatial_mean_mul(int coord_0, double *phi, double *psi){
  
  int i,j;
  int pos1=get_index(coord_0,0,0);
  int pos2=get_index(coord_0,0,0);
  double obsr=0.;
  
  for(i=0;i<Nx;i++){
    
    pos1=pos2;
    for(j=0;j<Ny;j++){
      
      obsr+=DObs(pos1)*psi[pos1];
      pos1=d_up[pos1*d+2];
    }
    pos2=d_up[pos2*d+1];
  }
  
  return obsr/((double)(Nx*Ny));
  
}





/* SPATIAL_MEAN_MUL: d-1 dimensional mean, used for TWI                                                                       */
/* There is no difference with respect the previously shown functions, with exception of the observables we chose to multiply */
/* In this case, this function is used to implement the chain rule of the forward lattice derivative for the rhs              */



double spatial_mean_mul2(int coord_0, double *phi){
  
  int i,j;
  int pos1=get_index(coord_0,0,0);
  int pos2=get_index(coord_0,0,0);
  double obsr=0.;
  
  for(i=0;i<Nx;i++){
    
    pos1=pos2;
    for(j=0;j<Ny;j++){
      
      obsr+=DObs(pos1)*phi[d_up[pos1*d]];
      pos1=d_up[pos1*d+2];
    }
    pos2=d_up[pos2*d+1];
  }
  
  return obsr/((double)(Nx*Ny));
  
}





/* SPATIAL_MEAN_MUL: d-1 dimensional mean, used for TWI                                                                       */
/* There is no difference with respect the previously shown functions, with exception of the observables we chose to multiply */
/* In this case, this function is used to implement the chain rule of the backward lattice derivative for the rhs             */



double spatial_mean_mul3(int coord_0, double *phi){
  
  int i,j;
  int pos1=get_index(coord_0,0,0);
  int pos2=get_index(coord_0,0,0);
  double obsr=0.;
  
  for(i=0;i<Nx;i++){
    
    pos1=pos2;
    for(j=0;j<Ny;j++){
      
      obsr+=DObs(pos1)*phi[d_down[pos1*d]];
      pos1=d_up[pos1*d+2];
    }
    pos2=d_up[pos2*d+1];
  }
  
  return obsr/((double)(Nx*Ny));
  
}




/* PHI2: define the mean of a squared lattice quantity */



double phi2(double *phi){

  int i;
  double sum=0.;

  for(i=0;i<Vol;i++) sum+=phi[i]*phi[i];

  return sum;

}





/* PHI4: define the mean of the fourth power of a lattice quantity */



double phi4(double *phi){

  int i;
  double sum=0.;

  for(i=0;i<Vol;i++) sum+=phi[i]*phi[i]*phi[i]*phi[i];

  return sum/(double)Vol;

}





/* PHI6: define the mean of a sixth power of a lattice quantity */



double phi6(double *phi){

  int i;
  double sum=0.;

  for(i=0;i<Vol;i++) sum+=phi[i]*phi[i]*phi[i]*phi[i]*phi[i]*phi[i];

  return sum;

}






/* PHI_TILDA: define the fourier transform of the field                                 */
/* L-for every lattice site                                                             */
/*  -take the value of the coordinate for the lattice index value                       */
/*  -calculate the real part of the field fourier transform as phi[i]*cos(p_j*x_j)      */
/*  -calculate the imaginary part of the field fourier transform as phi[i]*sin(p_j*x_j) */



void phi_tilda(double* p, double* phi,complex* phitilda){
  
  int i;
  static int x[3];
  phitilda->Re=phitilda->Im=0.0;
  
  for(i=0;i<Vol;i++){    
    
    coordinates(i,x);
    phitilda->Re+=(phi[i])*cos(p[0]*x[0]+p[1]*x[1]+p[2]*x[2]);
    phitilda->Im+=(phi[i])*sin(p[0]*x[0]+p[1]*x[1]+p[2]*x[2]);
    
  }
}





/* ACTION: define the action of the system                           */
/* L-for avery lattice site                                          */
/*  L-for every direction mu                                         */
/*   -evaluate the phi^4 local action                                */



double action(double *phi){
  
  double local=0.0;
  int mu,i;
  
  for(i=0;i<Vol;i++){
    
    local += (mass_square*0.5+d)*((Phi[i])*(Phi[i]))+lambda/24.0*(Phi[i]*Phi[i]*Phi[i]*Phi[i]);
    
    for(mu=0;mu<d;mu++)
   
   local += -(Phi[i])*(Phi[d_up[i*d+mu]]+Phi[d_down[i*d+mu]])*0.5;
 
  }
  
  return local/Vol;
  
}





/* EMT: function for the evalutation of the EMT components on lattice                                 */
/* symmetric formulation of lattice derivatives has been used: take a look to Scaracciolo and friends */



double EMT(int mu, int nu, int i, double *Field){
  
  if(mu!=nu){
    
    return (Field[d_up[i*d+mu]]-Field[d_down[i*d+mu]])*(Field[d_up[i*d+nu]]-Field[d_down[i*d+nu]])/4.0;
    
  }else{
    
    int rho;    
    double lagrangian=(mass_square+d*0.5)*0.5*((Field[i])*(Field[i]))+lambda/24.0*((Field[i])*(Field[i])*(Field[i])*(Field[i]));
    
    for(rho=0;rho<d;rho++){
      
      lagrangian +=  -0.25*(Field[d_up[i*d+rho]]*Field[d_down[i*d+rho]]);
      
    }
    
    return (Field[d_up[i*d+mu]]-Field[d_down[i*d+mu]])*(Field[d_up[i*d+nu]]-Field[d_down[i*d+nu]])/4.0-lagrangian;
    
  }
}





/* Zrhs: r.h.s of the integrated TWI                                                                                          */
/* This function evaluates the product of a given observable (defined in the macro above) and its time derivative             */
/* Both factors have to be evaluated in two different points along the time direction.                                        */
/* Being more precise, i0 must lie outside the integration domain with center in y0 and that characterizes the integrated TWI */
/* More precise informations can be found in the article of Del Debbio, Patella, Rago [http://arxiv.org/abs/1306.1173]        */


double Zrhs(int i0,int j0,double *phi,char name){
  
  if(name=='s'||name=='c')  return spatial_mean(i0,phi)*(spatial_mean(j0+1,phi)-spatial_mean(j0-1+Nt,phi))/2.0;
  else if(name=='f') return spatial_mean(i0,phi)*(spatial_mean(j0+1,phi)-spatial_mean(j0,phi));
  else if(name=='b') return spatial_mean(i0,phi)*(spatial_mean(j0,phi)-spatial_mean(j0-1+Nt,phi));
  else {
    
    mylog("[OBSERVABLES]: Need to specify a lattice derivative for Zrhs \n");
    exit(1);
    
  }  
}





/* VIOLATION: function that measure the numerical difference from naive and chain rule lattice derivative formulations*/



double violation(int i0, int j0, double *phi,char name){

  if(name=='s'||name=='c')   return -spatial_mean(i0,phi)*spatial_mean(j0+1,phi)/2.0+spatial_mean(i0,phi)*spatial_mean(j0-1+Nt,phi)/2.0+spatial_mean(i0,phi)*(spatial_mean_mul2(j0,phi)-spatial_mean_mul2(j0-1+Nt,phi))/2.0;
  else if(name=='f') return -spatial_mean(i0,phi)*spatial_mean(j0+1,phi)-spatial_mean(i0,phi)*spatial_mean(j0,phi)+spatial_mean(i0,phi)*spatial_mean_mul2(j0,phi);
  else if(name=='b') return -spatial_mean(i0,phi)*spatial_mean(j0,phi)-spatial_mean(i0,phi)*spatial_mean(j0-1+Nt,phi)+spatial_mean(i0,phi)*spatial_mean_mul2(j0-1+Nt,phi);
  else {
    
    mylog("[OBSERVABLES]: Need to specify a lattice derivative for violation check \n");
    exit(1);
    
  }  
}





/* ZRHS: r.h.s of the integrated TWI with chain rule lattice derivative */



double Zrhs_chain(int i0,int j0,double *phi,double *dphi){
  
  return spatial_mean(i0,phi)*spatial_mean_mul(j0,phi,dphi); 
  
}





/* Zlhs: l.h.s. of the integrated TWI                                                                                                                                                                     */
/* This function evaluates the product of a given observable (defined in the macro above) and the wall average of the variation of such observables times the integrated action of the Jacobian operator  */
/* Both factors have to be evaluated in two different points along the time direction                                                                                                                     */
/* Being more precise, i0 must lie outside the integration domain, with center in y0, that characterizes the integrated TWI                                                                               */
/* More precise informations can be found in the article of Del Debbio, Patella, Rago [http://arxiv.org/abs/1306.1173]                                                                                    */



double Zlhs(int i0, int j0, double *dphi, double *psi){

  return spatial_mean(i0,dphi)*spatial_mean_mul(j0,dphi,psi);


}





/* SPATIAL_MEAN: d-1 dimensional mean, used for TWI */
/* This function will give you the d-1 dimensional mean of a product of given observables (defined in the macro above), in a given point along the time direction */



double LHpert(int coord_0,int coord_1, double *phi, double *psi){
  
  int i,j;
  int pos1=get_index(coord_0,0,0);
  int pos2=get_index(coord_1,0,0);
  
  double obsr=0.;
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      
      obsr+=DObs(pos2)*psi[pos2]*Obs(pos1);
      pos1=d_up[pos1*d+2];
      pos2=d_up[pos2*d+2];
    }
    pos1=d_up[pos1*d+1];
    pos2=d_up[pos2*d+1];
  }
  
      return obsr/((double)(Nx*Ny));
  
}

double LHpert_exact(int coord_0,int coord_1,int size,int t, double *dphi, double *phi, double *jac_vec){
  
  double lh=0.;
  double temp=0.0;
  int i,j,k,l,m;
  int pos1,pos2;
  for(l=0;l<Nx;l++){
    for(m=0;m<Ny;m++){
      
      temp=0.;
      for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	  for(k=Nt/2-size;k<=Nt/2+size;k++){
	    pos1=get_index(abs(k-coord_1),abs(i-l),abs(j-m));
	    pos2=get_index(k,i,j);
	    temp+=jac_vec[pos1+Vol*t]*dphi[pos2];
	    
	  }
      temp*=Obs(get_index(coord_0,l,m))*DObs(get_index(coord_1,l,m));
      
      lh+=temp;
    }
  }
   
      return lh/((double)(Nx*Ny));
  
}


double LHpert_exact_vol(int dis,int size,int t, double *dphi, double *phi, double *jac_vec){
  
  double lh=0.;
  double temp=0.0;
  int i,j,k,l,m,n=0;
  int pos1,pos2;
  for(n=0;n<Nt;n++)
    for(l=0;l<Nx;l++){
      for(m=0;m<Ny;m++){
      
      temp=0.;
      for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	  for(k=0;k<=2*size;k++){
	    pos1=get_index(abs(k-size)%Nt,abs(i-l),abs(j-m));
	    pos2=get_index(abs(k+dis-size+n)%Nt,i,j);
	    temp+=jac_vec[pos1+Vol*t]*dphi[pos2];
	    
	  }
      temp*=Obs(get_index(n,l,m))*DObs(get_index((dis+n)%Nt,l,m));
      
      lh+=temp;
    }
  }
   
    return lh/(double)Vol;
  
}

