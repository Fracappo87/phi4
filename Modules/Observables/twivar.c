#include "twivar.h"
double delta_pr_1(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=dB_Phi(y,alpha)*Jacobian(tflow,DIST);

double sum2=0.;
sum2+=dB_Jacobian(tflow,DIST,alpha)*Phi(y);

return (sum1+sum2)*dB_Phiboundary(x,beta);
}

double delta_pr_2(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=sB_dC_dB_dB_Phi(y,alpha)*Jacobian(tflow,DIST);

double sum2=0.;
sum2+=sB_dC_dB_dB_Jacobian(tflow,DIST,alpha)*Phi(y);

return (sum1+sum2)*dB_Phiboundary(x,beta);
}

double delta_pr_3(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=sB_dB_dB_Phi(y)*dB_Jacobian(tflow,DIST,alpha);

double sum2=0.;
sum2+=sB_dB_dB_Jacobian(tflow,DIST)*dB_Phi(y,alpha);

return (sum1+sum2)*dB_Phiboundary(x,beta);
}

double delta_pr_4(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
int a11;

for(a11=0;a11<d;a11++){
sum1+=dC_dB_Phi(y,a11,alpha)*dB_Jacobian(tflow,DIST,a11);

}
double sum2=0.;
int a12;

for(a12=0;a12<d;a12++){
sum2+=dC_dB_Jacobian(tflow,DIST,a12,alpha)*dB_Phi(y,a12);

}
return (sum1+sum2)*dB_Phiboundary(x,beta);
}

double delta_pr_5(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=dB_Phi(y,alpha)*Phi(y)*Phi(y)*Jacobian(tflow,DIST);

sum1/=6.;

sum1*=3;

double sum2=0.;
sum2+=dB_Jacobian(tflow,DIST,alpha)*Phi(y)*Phi(y)*Phi(y);

sum2/=6.;

return (sum1+sum2)*dB_Phiboundary(x,beta);
}

double delta_pr_6(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=sB_dC_dB_dB_Phi(y,alpha)*Phi(y)*Phi(y)*Jacobian(tflow,DIST);

sum1/=6.;

sum1*=3;

double sum2=0.;
sum2+=sB_dC_dB_dB_Jacobian(tflow,DIST,alpha)*Phi(y)*Phi(y)*Phi(y);

sum2/=6.;

return (sum1+sum2)*dB_Phiboundary(x,beta);
}

double delta_pr_7(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=sB_dB_dB_Phi(y)*dB_Phi(y,alpha)*Phi(y)*Jacobian(tflow,DIST);

sum1/=6.;

sum1*=2;

double sum2=0.;
sum2+=sB_dB_dB_Phi(y)*dB_Jacobian(tflow,DIST,alpha)*Phi(y)*Phi(y);

sum2/=6.;

double sum3=0.;
sum3+=sB_dB_dB_Jacobian(tflow,DIST)*dB_Phi(y,alpha)*Phi(y)*Phi(y);

sum3/=6.;

return (sum1+sum2+sum3)*dB_Phiboundary(x,beta);
}

double delta_pr_8(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
int a11;

for(a11=0;a11<d;a11++){
sum1+=dC_dB_Phi(y,a11,alpha)*dB_Phi(y,a11)*Phi(y)*Jacobian(tflow,DIST);

}
sum1/=6.;

sum1*=2;

double sum2=0.;
int a12;

for(a12=0;a12<d;a12++){
sum2+=dC_dB_Phi(y,a12,alpha)*dB_Jacobian(tflow,DIST,a12)*Phi(y)*Phi(y);

}
sum2/=6.;

double sum3=0.;
int a13;

for(a13=0;a13<d;a13++){
sum3+=dC_dB_Jacobian(tflow,DIST,a13,alpha)*dB_Phi(y,a13)*Phi(y)*Phi(y);

}
sum3/=6.;

return (sum1+sum2+sum3)*dB_Phiboundary(x,beta);
}

double delta_pr_9(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
int a11;

for(a11=0;a11<d;a11++){
sum1+=dB_Phi(y,a11)*dB_Phi(y,a11)*dB_Phi(y,alpha)*Jacobian(tflow,DIST);

}
sum1/=6.;

double sum2=0.;
int a12;

for(a12=0;a12<d;a12++){
sum2+=dB_Phi(y,a12)*dB_Phi(y,a12)*dB_Jacobian(tflow,DIST,alpha)*Phi(y);

}
sum2/=6.;

double sum3=0.;
int a13;

for(a13=0;a13<d;a13++){
sum3+=dB_Phi(y,a13)*dB_Jacobian(tflow,DIST,a13)*dB_Phi(y,alpha)*Phi(y);

}
sum3/=6.;

sum3*=2;

return (sum1+sum2+sum3)*dB_Phiboundary(x,beta);
}

double delta_pr_10(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=dB_Phi(y,alpha)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Jacobian(tflow,DIST);

sum1/=120.;

sum1*=5;

double sum2=0.;
sum2+=dB_Jacobian(tflow,DIST,alpha)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Phi(y);

sum2/=120.;

return (sum1+sum2)*dB_Phiboundary(x,beta);
}

double delta_pr_11(int x,int y,int beta,int alpha,int tflow,double *Phi,double *Phiboundary){

double sum1=0.;
sum1+=dB_Phi(y,alpha)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Jacobian(tflow,DIST);

sum1/=5040.;

sum1*=7;

double sum2=0.;
sum2+=dB_Jacobian(tflow,DIST,alpha)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Phi(y)*Phi(y);

sum2/=5040.;

return (sum1+sum2)*dB_Phiboundary(x,beta);
}

