/* UPDATE LIBRARY: contains the implementation of Metropolis and Cluster algorithm */


#include <time.h>
#include <math.h>
#include "update.h"





/* GLOBAL COUNTERS      */
/* -acceptance counter  */
/* -call counter        */



long int acc=0;
long int call=0;
long int call_clu=0;
long int acc_clu=0;




/* PROBABILITY RATIO                                                                                               */
/* -define the difference  between intial and trial state energies, which is given by { -L(phi*[i])+L((phi)[i]) }  */
/* IMPORTANT: contribution from, say, the nth site comes from phi[i] and phi[i+mu]: you must keep track of this    */




double action_diff(double *phip_i,int i){
  
  int mu;
  double local;
  
  local = (mass_square*0.5+d)*((*phip_i)*(*phip_i)-(iPhi(i))*(iPhi(i)))+(g/6.0)*((*phip_i)*(*phip_i)*(*phip_i)-iPhi(i)*iPhi(i)*iPhi(i))+(lambda/24.0)*((*phip_i)*(*phip_i)*(*phip_i)*(*phip_i)-iPhi(i)*iPhi(i)*iPhi(i)*iPhi(i))+(kappa/720.0)*((*phip_i)*(*phip_i)*(*phip_i)*(*phip_i)*(*phip_i)*(*phip_i)-iPhi(i)*iPhi(i)*iPhi(i)*iPhi(i)*iPhi(i)*iPhi(i));
  
  for(mu=0;mu<d;mu++){
  
  local += -((*phip_i)-iPhi(i))*(iPhi(d_up[i*d+mu])+iPhi(d_down[i*d+mu]));
    
  }

  return local;  
  
}





/* METROPOLIS ALGORITHM                                                                                                                */
/* L-for every lattice site(indexed by i)                                                                                              */
/*  -update calls counter                                                                                                              */
/*  -initialize 2 random numbers to be used for metropolis test                                                                        */
/*  -choose a new configuration as the new probe state                                                                                 */
/*  IF-exp(-dE) (probability ratio) > 1 the new state is accepted, acceptance counter updated                                          */
/*  IF-exp(-dE) (probability ratio) <  do metropolis test: if ran[1] < exp(-dE),  new state is accepted, acceptance counter updated,   */
/*  otherwise it is refused                                                                                                            */



void metropolis(double delta){

  int i;
  double phip_i;
  static double rand[2];
  
  for(i=0;i<Vol;i++){
    
    call++;
    ranlxd(rand,2);                                                       
    phip_i = iPhi(i)+(1-rand[0]*2.0)*delta;                                                   
    double dE=action_diff(&phip_i,i);                  
    
    if(rand[1]<exp(-dE)){                                                
      
      iPhi(i)=phip_i; 
      acc++;
    
    }
  }
}





/* CLUSTER VECTOR INITIALIZATION                                    */
/* allocate the memory for cluster and cluster weight vectors       */
/* consistency check on memory allocation                           */

static double *rand_flip;

void cluster_init(){
  
  cluster=(int*)malloc(Vol*sizeof(int));
  cluster_weight=(int*)malloc(Vol*sizeof(int)); 
  rand_flip=(double*)malloc(Vol*sizeof(double));

  
  if(cluster==NULL){
    
    mylog("[UPDATE]: Cluster vectors allocation error\n");
    mylog("[UPDATE]: Run stopped\n");
    exit(1);

  } 
  
  if(cluster_weight==NULL){
    
    mylog("[UPDATE]: Cluster_weight vector  allocation error\n");
    mylog("[UPDATE]: Run stopped\n");
    exit(1);
    
  }
  
}





/* CLUSTER PROBABILITY                                                                                    */
/* -Evaluates the probability of link creation between lattice sites, which is an ising-type probability  */



double cluster_prob(int i,int mu){

  return 1-exp(-imod_Phi(i)*imod_Phi(d_up[i*d+mu])-iPhi(i)*iPhi(d_up[i*d+mu]));

}





/* CLUSTER FIND ROOT                                                                                                                         */ 
/* W-until cluster[i] is different from i                                                                                                    */
/*  -set i_tmp=cluster[i]                                                                                                                    */
/*  IF-cluster[i] different from cluster[i_tmp] (that is i_tmp is not the root), set cluster[i]=cluster[i_tmp] (that is continue to search)  */
/*  -set i=i_tmp                                                                                                                             */



int find_root(int i){
  
  int i_tmp,j_tmp=i;
  
    while(cluster[j_tmp]!=j_tmp){
    
    i_tmp=cluster[j_tmp];
   
    if(cluster[j_tmp]!=cluster[i_tmp]){
      
      cluster[j_tmp]=cluster[i_tmp];
      
    }
    
    j_tmp=i_tmp;
    
    }
 
  return j_tmp;
}





/*CLUSTER ALGORITHM                                                                                                                                                                                                    */
/* L-for every lattice site                                                                                                                                                                                            */
/*  -set the initial cluster vector entries giving "-1" for every lattice site                                                                                                                                         */
/*  -set the intiali cluster_weight vector entries to 0                                                                                                                                                                */
/*  -initialize D-dimensional vector with random entries to be used for linking test                                                                                                                                   */
/* L-for avery lattice site                                                                                                                                                                                            */
/*  L-for every direction mu                                                                                                                                                                                           */
/*  IF-cluster probability > rand_cluster[mu], then create a link between the sites                                                                                                                                    */
/*    IF-both sites not linked, give them a root (i) to both sites, set cluster weight of (root) to 1                                                                                                                  */
/*    IF-one site is linked, the other not, then give the same root of the already linked site, update the cluster weight of (root)                                                                                    */
/*    IF-same of (3), reversed case, do the same (but in reverse way)                                                                                                                                                  */
/*    IF-both sites are linked                                                                                                                                                                                         */
/*      IF-cluster_weight[find_root(i)]<cluster_weight[find_root(i+mu)], then convert the first cluster, giving it the root of (i+mu), update the cluster weight(this operation must be done before changing the root) */
/*      ELSE-cluster_weight[find_root(i)]≥cluster_weight[find_root(i+mu)]                                                                                                                                              */
/*          IF cluster[i]!=cluster[i+mu], then convert the second cluster, giving it the root of (i), update the cluster weight(this operation must be done before changing the root)                                  */       
/*          ELSE-cluster[i]=cluster[i+mu], then we are linking a single cluster with itself, no root changing are needed, update cluster weight adding 1                                                               */



void cluster_creation(){

  int i,mu;
  double rand_cluster[d];

  for(i=0;i<Vol;i++){
  
  cluster[i]=-1;
  cluster_weight[i]=0;
  
  }   
  
  for(i=0;i<Vol;i++){
    
    ranlxd(rand_cluster,d);     
    
    for(mu=0;mu<d;mu++){
    
      call_clu++;
      
      if(rand_cluster[mu]<cluster_prob(i,mu)){   
	
	acc_clu++;
	
	if(cluster[i]==-1 && cluster[d_up[i*d+mu]]==-1){           
	  
	  cluster[i]=i;
	  cluster[d_up[i*d+mu]]=i;
	  cluster_weight[i]=1;
	  
	}
	else if(cluster[i]!=-1 && cluster[d_up[i*d+mu]]==-1){        
	  
	  cluster[i]=find_root(i);
	  cluster[d_up[i*d+mu]]=cluster[i];
	  cluster_weight[cluster[i]]+=1;
	  
	}
	else if(cluster[i]==-1 && cluster[d_up[i*d+mu]]!=-1){        
	
	  cluster[d_up[i*d+mu]]=find_root(d_up[i*d+mu]); 
	  cluster[i]=cluster[d_up[i*d+mu]];
	  cluster_weight[cluster[d_up[i*d+mu]]]+=1;
	  
	}
	else if(cluster[i]!=-1 && cluster[d_up[i*d+mu]]!=-1){                              
	  if(cluster_weight[find_root(i)]<cluster_weight[find_root(d_up[i*d+mu])]){  
	  	  
	    cluster_weight[find_root(d_up[i*d+mu])] += cluster_weight[find_root(i)]+1;   	    
	    cluster_weight[find_root(i)]=0;
	    cluster[find_root(i)]=find_root(d_up[i*d+mu]);
	    cluster[i]=cluster[find_root(i)];
	    
	  }else{
	    if(cluster[i]!=cluster[d_up[i*d+mu]]){
	    
	    cluster_weight[find_root(i)] += cluster_weight[find_root(d_up[i*d+mu])]+1;     	   
	    cluster_weight[find_root(d_up[i*d+mu])]=0;
	    cluster[find_root(d_up[i*d+mu])]=find_root(i);
	    cluster[d_up[i*d+mu]]=cluster[find_root(d_up[i*d+mu])];
	    
	    }else{

	      cluster_weight[find_root(i)] += 1;     	   
	      cluster[d_up[i*d+mu]]=find_root(i);
	      
	    }
	  }
	}
      }
    }
  }
}





/* CLUSTER FLIP: function implementing cluster sign flip                                      */
/* L-for every lattice site                                                                   */
/*  IF-cluster weight[i]>0, then we are considering a root with a non zero dimension cluster  */
/*    -choose a random number betwenn 0 and 1                                                 */
/*    IF-rand>0.5, then we decide to flip the cluster sign                                    */
/*      L-for every lattice site                                                              */
/*       IF-cluster[j]>=0,then we are choosign swtiched-on site, otherwise find_root crashes  */
/*         IF-find_root(i)=j, then the site cluster variables is the same of the cluster root */
/*           -flip the sign of the filed at site j, i.e. change the cluster fields' sign      */  
 

void cluster_flip(){
  int i;

  for(i=0; i<Vol;i++)
    if(cluster_weight[i]>0)
      ranlxd(rand_flip+i,1);
      
  for(i=0;i<Vol;i++){
      if(cluster[i]>=0)
	if(rand_flip[find_root(i)]<0.5)
	  iPhi(i)*=-1.0;
  }
  
}


/*
void cluster_flip(){

  int i,j;
  double rand_flip[1];
  
  for(i=0;i<Vol;i++){
    if(cluster_weight[i]>0){
      
      ranlxd(rand_flip,1);
      
      if(rand_flip[0]<0.5){
	for(j=0;j<Vol;j++){
	  if(cluster[j]>=0){
	   if(find_root(j)==i){
	    
	    iPhi(j)*=-1.0;
	    
	   }
	  }
	}
      }
    }
  }
}

*/



/* UPDATE_SWEEP: function implementing the sweep for our update algorithm */
/* L-for each sweep step                                                  */
/*  -create clusters                                                      */
/*  -flip clusters signs                                                  */
/*  L-for each sweep_ratio step                                           */
/*   -make a Metropolis step                                              */



void update_sweep(int N, int n){

  int i,l;
  

    for(l=0;l<N;l++){
     
      cluster_creation();   
      cluster_flip();

      for(i=0;i<n;i++) metropolis(delta);

    }  
 
}





/* CLUSTER FREE */
/* -deallocate cluster and cluster_weight vector memories */



void free_cluster(){

  free(rand_flip);
  free(cluster);
  free(cluster_weight);
  cluster=NULL;
  cluster_weight=NULL;

  mylog("[UPDATE]: cluster vector's memory freed\n");
  mylog("[UPDATE]: cluster_weight vector's memory freed\n");
  mylog("\n");

}

   
 
  
