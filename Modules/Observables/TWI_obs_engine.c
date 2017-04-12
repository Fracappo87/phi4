 /* TWI_OBS_ENGINE: routine for the definitions of all the operators that enters in the definition of our local translation Ward Identities*/



#include "TWI_obs_engine.h"
#include "probevec"





/* TWI_operator: vector containing all the operators that mix with the energy momentum tensor */



obs_ptr_t TWI_operator[NUM_OP]={phi_1,phi_2,phi_4,phi_5};





double *avg_vector_init(void){
  
  double *vec=(double*)malloc((NUM_OP+1)*NUM_PROB*(N_flowtimes+1)*sizeof(double));
  
  return vec;
  
}





void vec_clearer(double *vec){

  int Ntot=(NUM_OP+1)*NUM_PROB*(N_flowtimes+1);
  int j;

  for(j=0;j<Ntot;j++)
    vec[j]=0.;

}




/*     Energy momentum tensor mixing operators       */
/*                                                   */
/*Definitions   (mu,nu index not summed)             */
/*                                                   */
/*phi_1 = (d_muPhi*d_nuPhi)/2!                       */
/*                                                   */
/*phi_2 = Phi^2/2!                                   */
/*                                                   */
/*phi_3 = delta_{munu}(d_muPhi*d_muPhi)/2!           */
/*                                                   */
/*phi_4 = delta_{munu}(Phi^4)/4!                     */
/*                                                   */
/*phi_5 = delta_{munu}sum_rho[d_rhoPhi*d_rhoPhi]/2!  */
/*                                                   */
/*phi_6 = delta_{munu}(Phi^6)/6!                     */
/*                                                   */
/*phi_7 = delta_{munu}(PhiBOXPhi)/2!                 */
/*                                                   */
/*phi_8 = delta_{munu}(Phid_mud_muPhi)/2!            */



double phi_1(int i, int mu, int nu, double *Phi){

  return  dB_Phi(i,mu)*dB_Phi(i,nu)*0.5;  

}





double phi_2(int i, int mu, int nu, double *Phi){

  if(mu==nu)   return Phi[i]*Phi[i]/2.0;
  else return 0.;

}





double phi_3(int i, int mu, int nu, double *Phi){

  if(mu==nu) {

    double temp=dB_Phi(i,mu);
    return temp*temp/2.0;

  }
    else return 0.;

}





double phi_4(int i, int mu, int nu, double *Phi){

  if(mu==nu)  {
    double temp=Phi[i]*Phi[i];
    return temp*temp/24.0;
  }
    else return 0.;
  
}





double phi_5(int i, int mu, int nu, double *Phi){

  if(mu==nu){
    
    int rho;
    double sum=0.;
  
    for(rho=0;rho<d;rho++){

      double temp=dB_Phi(i,rho);
      sum+=temp*temp/2.0;

	}
    return sum;
  }else{

    return 0.;
  }
}





double phi_6(int i, int mu, int nu, double *Phi){

  if(mu==nu){
    double temp=Phi[i]*Phi[i]*Phi[i];
    return temp*temp/720.;
  }
    else return 0.;
}





double phi_7(int i, int mu, int nu, double *Phi){

  if(mu==nu){
    
    //int rho;
    //double sum=0.;
    
    //for(rho=0;rho<d;rho++)
    //sum+=dB_dB_Phi(i,rho);
    
    return sB_dB_dB_Phi(i)*Phi[i]/2.;//sum*Phi[i]/2.;
    
  }else{
    
    return 0.;
    
  }
}





double phi_8(int i, int mu, int nu, double *Phi){

  if(mu==nu){
    
    return Phi[i]*dB_dB_Phi(i,mu)/2.;
 
 }else{

    return 0.;

  }
}






double RHS_loc(int x,int y,int nprob,int ntens,double *phi_flow,double *phi){

  double res=0.;
  int mu,nu;
  double deltapr;
    for(mu=0;mu<d;mu++)
      for(nu=0;nu<d;nu++){
	
	
	//deltapr=(PROB_OP[nprob](d_up[y*d+mu],nu,phi_flow)-PROB_OP[nprob](d_down[y*d+mu],nu,phi_flow))/2.;
	
	deltapr=-(TWI_operator[ntens](d_up[x*d+mu],mu,nu,phi)-TWI_operator[ntens](d_down[x*d+mu],mu,nu,phi))/2.;
	
	//res+=TWI_operator[ntens](x,mu,nu,phi)*deltapr;                                   
	res+=PROB_OP[nprob](y,nu,phi_flow)*deltapr;                                         
	
}


  return res;
}





/* RHS: functions for the evaluation of the RHS of our TWI, taking rigid volume averages                                    */
/* It simply glues up all the pieces that are included in the definition of the RHS, summing correctly over Lorentz indices */
/* For more details about the implementation, look at Dropbox/Scalaer_EMT/Notes_Plymouth/TWI_RC.pdf                         */



double RHS(int x,int y,int nprob,int ntens,double *phi_flow,double *phi){

  int i,j,k;
  int pos=x,pos1=y;

  double corr=0.;
  for(i=0;i<Nt;i++){
    for(j=0;j<Nx;j++){
      for(k=0;k<Ny;k++){
	
        corr+=RHS_loc(pos,pos1,nprob,ntens,phi_flow,phi);
        pos=d_up[pos*d+2];
        pos1=d_up[pos1*d+2];
      }
      pos=d_up[pos*d+1];
      pos1=d_up[pos1*d+1];
    }
    pos=d_up[pos*d];
    pos1=d_up[pos1*d];
  }
  return corr/(double)(Vol);

}


/* LHS: functions for the evaluation of the RHS of our TWI, taking rigid volume averages                                    */
/* It simply glues up all the pieces that are included in the definition of the RHS, summing correctly over Lorentz indices */
/* For more details about the implementation, look at Dropbox/Scalaer_EMT/Notes_Plymouth/TWI_RC.pdf                         */



double LHS_delta_V(int nprob,int x,int y, int tflow, double *phi_flow,double *phi){

  int i,j,k,nu;
  int pos=x,pos1=y;
  
  double corr=0.;
  for(i=0;i<Nt;i++){
    for(j=0;j<Nx;j++){
      for(k=0;k<Ny;k++){
	for(nu=0;nu<d;nu++)
	  corr+=Delta_PROB_OP[nprob](pos,pos1,nu,nu,tflow,phi_flow,phi);
	pos=d_up[pos*d+2];
        pos1=d_up[pos1*d+2];
      }
      pos=d_up[pos*d+1];
      pos1=d_up[pos1*d+1];
    }
    pos=d_up[pos*d];
    pos1=d_up[pos1*d];
  }
  return corr/(double)(Vol);

}








