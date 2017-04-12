/* GLOB_CONST: library for global constant definitions*/



#ifndef __GLOB_CONST
#define __GLOB_CONST





/* MACRO: used to define global constants, vectors and structures, if main program is defined */



#ifdef MAIN_PROGRAM
#  define GLB_VAR(type,name,init...) type name init
#else
#  define GLB_VAR(type,name,init...) extern type name
#endif





/* GREEK PI, square mas... */



#define mPI 3.1415926535897932385
//#define mass_square -mass*mass





/* GLOBAL VARIABLES*/



GLB_VAR(int,Nt,=0);                     /* local lattice size in direction T       */ 
GLB_VAR(int,Nx,=0);                     /* local lattice size in direction X       */
GLB_VAR(int,Ny,=0);                     /* local lattice size in direction Y       */
GLB_VAR(int,Nz,=0);                     /* local lattice size in direction Z       */
GLB_VAR(int,Vol,=0);                    /* volume                                  */
GLB_VAR(int,d,=0);                      /* dimension of the system Z               */
GLB_VAR(int,d_field,=0);                /* number of field components              */
GLB_VAR(int,N_traj,=0);                 /* number of trajectories                  */
GLB_VAR(int,N_sweeps,=0);               /* number of measures                      */
GLB_VAR(int,sweep_ratio,=0);            /* number of measures                      */
GLB_VAR(int,N_term,=0);                 /* number of term steps                    */
GLB_VAR(int,N_report,=0);               /* report rate                             */
GLB_VAR(int,N_save,=0);                 /* save rate                               */
GLB_VAR(int,Last_count,=0);             /* number of the initial configuration     */
GLB_VAR(int,N_meas,=0);                 /* measure rate                            */
GLB_VAR(double,delta,=0);               /* metropolis controlo parameter           */
GLB_VAR(double,mass,=0);                /* square mass parameter                   */
GLB_VAR(double,mass_square,=0);         /* mass parameter                          */
GLB_VAR(double,lambda,=0);              /* phi^4 coupling                          */
GLB_VAR(double,g,=0);                   /* phi^3 coupling                          */
GLB_VAR(double,kappa,=0);               /* phi^6 coupling                          */
GLB_VAR(int,int_ord,=0);                /* integrator order                        */
GLB_VAR(double,Bulk_mass,=0);           /* Bulk mass parameter                     */
GLB_VAR(double,Bulk_lambda,=0);         /* Bulk phi^4 coupling                     */
GLB_VAR(double,Bulk_g,=0);              /* Bulk phi^3 coupling                     */
GLB_VAR(double,Bulk_kappa,=0);          /* Bulk phi^6 coupling                     */
GLB_VAR(int,seed,=0);                   /* randomseed                              */
GLB_VAR(int *,d_up,=NULL);              /* nnbr up vector                          */
GLB_VAR(int *,d_down,=NULL);            /* nnbr down vector                        */
GLB_VAR(double *,Phi,=NULL);            /* field vector                            */
GLB_VAR(int *,cluster,=NULL);           /* cluster vector                          */
GLB_VAR(int *,cluster_weight,=NULL);    /* cluster wight vector                    */
GLB_VAR(double,dt,=0);                  /* flow time integration interval          */
GLB_VAR(int,N_steps,=0);                /* number of integration steps             */
GLB_VAR(int,N_flowtimes,=0);            /* number of flow times for measuring      */
GLB_VAR(int,N_flow,=0);                 /* measure rate during flow evolution      */
GLB_VAR(int *,Flow_vec,=NULL);          /* flow times vector                       */
GLB_VAR(char*,status,=NULL);            /* status variable (hot, cold, file start) */
GLB_VAR(char*,mass2_sign,=NULL);        /* square mass sign variable               */
GLB_VAR(int,traslpos,=0);               /* traslation point index                  */
GLB_VAR(int,probepos,=0);               /* probe point index                       */
GLB_VAR(int,rhs_measure,=0);            /* do you want to measure the rhs ?  (1,0) */
GLB_VAR(int,lhs_measure,=0);            /* do you want to measure the lhs ?  (1,0) */
GLB_VAR(int,avg_TWI,=0);                /* block lenght for intermediate averages  */
GLB_VAR(int,ndistances,=0);             /*Number of distances for TWI when MULTIDIS is defined */
#undef GLB_VAR






#endif

/* typedef struct field{ */

/*   double comp; */

/* }Field; */



