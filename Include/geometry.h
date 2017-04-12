/* GEOMETRY: library for system geometry initialization */



#ifndef __GEOMETRY
#define __GEOMETRY

#include"glob_const.h"
#include"input_output.h"
#include<stdlib.h>





/* COORDINATES */



void coordinates(int i, int* x);





/* GET_INDEX */



int get_index(int t,int x, int y);





/* COORDVEC */
typedef struct coord { int x[3]; } coord;
coord *coordvec;





/* GEOMETRY ALGORITHM  */



void  geometry();





/* FREE_GEO */



void free_geo();





#endif
