/* INPUT_OUTPUT: library for input/output function definition and implementation*/



#ifndef __INPUT_OUTPUT
#define __INPUT_OUTPUT


#include"glob_const.h"
#include<string.h> 
#include<stdarg.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/types.h>
#include<dirent.h>
#include<stdbool.h>





FILE* stream;



/* INIT_MYLOG */



void init_mylog(char *my_FILE, int narg);





/* MYLOG */



void mylog( const char * report_info, ... );





/* SEARCH_STRING */



int searchword(FILE* SourceFile,char* searchingString);





/* ERROR */



void error(bool stat, char* statString);



/* READINPUT */



void readinput(); 





/* FLOWINPUT */



void flow_input(char *my_FILE); 





/* TWI_INPUT */



void TWI_input(char*my_FILE);
void TDWI_input(char*my_FILE);





/* WRITE_OUTPUT */



void write_output(char *my_FILE, double *Obs, double Obs_check);





/* SEED_READINPUT */



void seed_readinput(char *my_FILE); 





/* SEED_OUTPUT */



void seed_output(char *my_FILE);








#endif
