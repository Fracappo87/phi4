#include "twiprobe.h"
double probe_1(int x, int alpha,double *Phi) {

double sum=0.;
sum+=dB_Phi(x,alpha)*Phi(x);

return sum;

}

double probe_2(int x, int alpha,double *Phi) {

double sum=0.;
sum+=sB_dC_dB_dB_Phi(x,alpha)*Phi(x);

return sum;

}

double probe_3(int x, int alpha,double *Phi) {

double sum=0.;
sum+=sB_dB_dB_Phi(x)*dB_Phi(x,alpha);

return sum;

}

double probe_4(int x, int alpha,double *Phi) {

double sum=0.;
int a1;

for(a1=0;a1<d;a1++){
sum+=dC_dB_Phi(x,a1,alpha)*dB_Phi(x,a1);

}
return sum;

}

double probe_5(int x, int alpha,double *Phi) {

double sum=0.;
sum+=dB_Phi(x,alpha)*Phi(x)*Phi(x)*Phi(x);

return sum/6.;

}

double probe_6(int x, int alpha,double *Phi) {

double sum=0.;
sum+=sB_dC_dB_dB_Phi(x,alpha)*Phi(x)*Phi(x)*Phi(x);

return sum/6.;

}

double probe_7(int x, int alpha,double *Phi) {

double sum=0.;
sum+=sB_dB_dB_Phi(x)*dB_Phi(x,alpha)*Phi(x)*Phi(x);

return sum/6.;

}

double probe_8(int x, int alpha,double *Phi) {

double sum=0.;
int a1;

for(a1=0;a1<d;a1++){
sum+=dC_dB_Phi(x,a1,alpha)*dB_Phi(x,a1)*Phi(x)*Phi(x);

}
return sum/6.;

}

double probe_9(int x, int alpha,double *Phi) {

double sum=0.;
int a1;

for(a1=0;a1<d;a1++){
sum+=dB_Phi(x,a1)*dB_Phi(x,a1)*dB_Phi(x,alpha)*Phi(x);

}
return sum/6.;

}

double probe_10(int x, int alpha,double *Phi) {

double sum=0.;
sum+=dB_Phi(x,alpha)*Phi(x)*Phi(x)*Phi(x)*Phi(x)*Phi(x);

return sum/120.;

}

double probe_11(int x, int alpha,double *Phi) {

double sum=0.;
sum+=dB_Phi(x,alpha)*Phi(x)*Phi(x)*Phi(x)*Phi(x)*Phi(x)*Phi(x)*Phi(x);

return sum/5040.;

}

