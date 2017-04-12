/* UPDATE LIBRARY: contains the implementation of Metropolis and Cluster algorithm */



#ifndef __UPDATE
#define __UPDATE


#include<math.h>
#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include"glob_const.h"
#include"lattice_config.h"
#include"geometry.h"
#include"ranlux.h"
#include"input_output.h"




/* GLOBAL COUNTERS      */
/* -acceptance counter  */
/* -call counter        */



long int acc,call,call_clu,acc_clu;





/* ACTION_DIFF */



double action_diff(double * phip_i,int i);

                         



/* METROPOLIS ALGORITHM */



void metropolis(double delta);





/* CLUSTER VECTORS INITIALIZATION */



void cluster_init();





/* CLUSTER PROBABILITY */



double cluster_prob();





/* CLUSTER FIND ROOT */



int find_root(int i);





/* CLUSTER ALGORITHM */



void cluster_creation();





/* CLUSTER FLIP */



void cluster_flip();
void cluster_flip2();




/* UPDATE_SWEEP */



void update_sweep(int N, int n);





/* CLUSTER FREE */



void free_cluster();





#endif
