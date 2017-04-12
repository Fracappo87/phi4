 /* TWI_OBS_ENGINE: routine for the definitions of all the operators that enters in the definition of our local translation Ward Identities*/
#include <math.h>
#include <stdio.h>
#include "DWI_obs_engine.h"
#include "observables.h"
#include "geometry.h"
#include "lattice_config.h"
#include "Jacobian.h"

/* TWI_operator: vector containing all the operators that mix with the energy momentum tensor */

dwiobs_ptr_t DWI_operator[dwiNUM_OP]={dwiphi_2,dwiphi_4,dwiphi_5,dwiphi_6,dwiphi_0};

dwiprb_ptr_t dwiPROB_OP[dwiNUM_PROB]={dwipr_1,dwipr_2,dwipr_3,dwipr_4,dwipr_5,dwipr_6,dwipr_7,dwipr_8};

dwiprb_ptr_t DPHI_PROB_OP[dwiNUM_PROB]={dphi_pr_1,dphi_pr_2,dphi_pr_3,dphi_pr_4,dphi_pr_5,dphi_pr_6,dphi_pr_7,dphi_pr_8};

dwiprb_der_ptr_t DDMUPHI_PROB_OP[dwiNUM_PROB]={ddmuphi_pr_1,ddmuphi_pr_2,ddmuphi_pr_3,ddmuphi_pr_4,ddmuphi_pr_5,ddmuphi_pr_6,ddmuphi_pr_7,ddmuphi_pr_8};

/* Norm_factor: vector containing all numerical coefficients needed to define the variations of our flowed probes */

double dim_factor[dwiNUM_PROB]={1.,2.,3.,4.,5.,0.,3.,4.};

double *dwiavg_vector_init(void){
  
  double *vec=malloc((dwiNUM_OP+1)*dwiNUM_PROB*(N_flowtimes+1)*sizeof(double));
  
  return vec;
  
}


void dwivec_clearer(double *vec){

  int Ntot=(dwiNUM_OP+1)*dwiNUM_PROB*(N_flowtimes+1);
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

double dwiphi_2(int i, double *phi){

  return phi[i]*phi[i]/2.0;
  
}


double dwiphi_4(int i,  double *phi){


    double temp=phi[i]*phi[i];
    return temp*temp/24.0;
  
}


double dwiphi_5(int i, double *Phi){


    
    int rho;
    double sum=0.;
  
    for(rho=0;rho<d;rho++){
      double temp=dB_Phi(i,rho);
      sum+=temp*temp;
	}
    return sum/2.;

}


double dwiphi_6(int i, double *phi){

 
    double temp=phi[i]*phi[i]*phi[i];
    return temp*temp/720.;
 
}


double dwiphi_0(int i, double *phi){

    return 1.;
 
}

/*   Probe observables definitions                */
/*   The basic construct is given by:             */
/*                                                */
/*   V_mu=1/n! d_mu phi(t,y)^n =                  */
/*       =1/(n-1)! phi(t,y)^(n-1)d_mu phi(t,y)    */
/*                                                */
/*  For more details, read the notes in:          */
/*  /Dropbox/Scalar_EMT/Notes_Plymouth/TWI_RC.pdf */


double dwipr_1(int i,double *phi){
  return phi[i]*phi[i]/2.;

}


double dwipr_2(int i,double *phi){
  
  double phi2=phi[i]*phi[i];
  return phi2*phi2/24.;
 
}


double dwipr_3(int i,double *phi){

  double phi3=phi[i]*phi[i]*phi[i];
  return phi3*phi3/720.;

}


double dwipr_4(int i,double *phi){

  double phi2=phi[i]*phi[i];
  double phi4=phi2*phi2;

  return phi4*phi4/40320.;
}


double dwipr_5(int i,double *phi){

  double phi2=phi[i]*phi[i];
  double phi4=phi2*phi2;

  return phi4*phi4*phi2/3628800.;
}


double dwipr_6(int i,double *phi){

  return 1.;

}


double dwipr_7(int i, double *Phi){

  int rho;
  double sum=0.;
  
  for(rho=0;rho<d;rho++){
    double temp=dB_Phi(i,rho);
    sum+=temp*temp;
  }
  return sum/2.;

}


double dwipr_8(int i,double *Phi){

  int rho;
  double sum=0.;
  
  for(rho=0;rho<d;rho++){
    double temp=dB_Phi(i,rho);
    sum+=temp*temp;
  }
  return sum*Phi[i]*Phi[i]/4.;

}


double dphi_pr_1(int i,double *phi){
  return phi[i];

}


double dphi_pr_2(int i,double *phi){
  
  double phi2=phi[i]*phi[i];
  return phi2*phi[i]/6.;
 
}


double dphi_pr_3(int i,double *phi){

  double phi3=phi[i]*phi[i]*phi[i];
  return phi3*phi[i]*phi[i]/120.;

}


double dphi_pr_4(int i,double *phi){

  double phi2=phi[i]*phi[i];
  double phi4=phi2*phi2;

  return phi4*phi2*phi[i]/5040.;
}


double dphi_pr_5(int i,double *phi){

  double phi2=phi[i]*phi[i];
  double phi4=phi2*phi2;

  return phi4*phi4*phi[i]/362880.;
}


double dphi_pr_6(int i,double *phi){

  return 0.;

}


double dphi_pr_7(int i,double *phi){

  return 0.;

}


double dphi_pr_8(int i,double *Phi){

int rho;
  double sum=0.;
  
  for(rho=0;rho<d;rho++){
    double temp=dB_Phi(i,rho);
    sum+=temp*temp;
  }
  return sum*Phi[i]/2.;


}



double ddmuphi_pr_1(int i,int rho,double *phi){
  return 0.;
}


double ddmuphi_pr_2(int i,int rho,double *phi){
  return 0.;
}


double ddmuphi_pr_3(int i,int rho,double *phi){
  return 0.;
}


double ddmuphi_pr_4(int i,int rho,double *phi){
  return 0.;
}


double ddmuphi_pr_5(int i,int rho,double *phi){
  return 0.;
}


double ddmuphi_pr_6(int i,int rho,double *phi){
  return 0.;
}


double ddmuphi_pr_7(int i,int rho,double *Phi){
  return dB_Phi(i,rho);
}


double ddmuphi_pr_8(int i,int rho,double *Phi){
  return dB_Phi(i,rho)*Phi[i]*Phi[i]/2.;
}


double dt_phi(int i, int tflow, double *phi){
  
  int x[d],y[d],xy[d];
  coordinates(i,x);
  
  
  int j,jp;
  double der=0.;
  for(j=0;j<Vol;j++){
    coordinates(j,y);
    for(jp=0;jp<d;jp++) xy[jp]=abs(x[jp]-y[jp]);
    
    der+=dt_jac_v[get_index(xy[0],xy[1],xy[2])+Vol*tflow]*phi[j];
  }
  return der;
}


double dtdmu_phi(int i,int mu,int tflow,double *Phi){
  int x[d],y[d],xy[d];
  coordinates(i,x);
  
  
  int j,jp;
  double der=0.;
  for(j=0;j<Vol;j++){
    coordinates(j,y);
    for(jp=0;jp<d;jp++) xy[jp]=abs(x[jp]-y[jp]);
   
    der+=dt_jac_v[get_index(xy[0],xy[1],xy[2])+Vol*tflow]*dB_Phi(j,mu);;
  }
  return der;




}
/*DELTA_V: function that evaluates the local infinitesimal translation of the probe.                     */
/* In order to make the evaluation as fastes as possible, we used a tabulation of the Jacobian operator, */
/* as weel as of the vectors that identify the points of probe position and translation application      */
/* For more details about the implementation, look at Dropbox/Scalaer_EMT/Notes_Plymouth/TWI_RC.pdf      */


double dwiLHS(int nprob,int y, int tflow, double *phi_flow, double *phi){
  
  double term1=0.;
  double term2=dim_factor[nprob]*dwiPROB_OP[nprob](y,phi_flow);
  double term3=0.;

  if(tflow==0){
    term1=0.;
    term3=0.;
  }else{
    
      term1=2.*(double)(tflow)*dt*DPHI_PROB_OP[nprob](y,phi_flow)*dt_phi(y,tflow,phi);
      int mu;
      for(mu=0;mu<d;mu++) term3+=2.*(double)(tflow)*dt*DDMUPHI_PROB_OP[nprob](y,mu,phi_flow)*dtdmu_phi(y,mu,tflow,phi);

  }


  return term1+term2+term3;


}


double dwiRHS(int nprob,int ntens,int y,double *phi_flow,double *phi){

  int i;
  double temp=0.;

  for(i=0;i<Vol;i++){

    temp+=DWI_operator[ntens](i,phi);

  }

  return dwiPROB_OP[nprob](y,phi_flow)*temp;

}
