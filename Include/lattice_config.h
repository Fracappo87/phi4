/* LATTICE_CONFIG: configures the initial values of fields in every lattice site */



#ifndef __LATTICE_CONFIG
#define __LATTICE_CONFIG

#include"glob_const.h"
#include<stdlib.h>





/* MACRO */ 
/* -define iPhi(site) as the component of the field at each lattice site          */
/* -define imod_Phi(site) as the abs of the field component at each lattice site  */



//#define iPhi(site) Phi[(site)].comp
#define iPhi(site) Phi[(site)]
//#define imod_Phi(site) sqrt(Phi[(site)].comp*Phi[(site)].comp)
#define imod_Phi(site) sqrt(Phi[(site)]*Phi[(site)])





/* CONFIGURATION */



void configuration();





/* FREE_LATTICE */




void free_Lattice();









#endif
