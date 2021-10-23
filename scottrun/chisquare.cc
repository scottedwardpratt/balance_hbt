#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
//#include <conio.h>
#include "parametermap.h"
#include <stdlib.h>


double  CalcChiSquared_pion(double &A){
	FILE *fptr=fopen("results/pipluspiplus/cf_outsidelong.dat","r");
	FILE *fptr1=fopen("cf_exp.dat","r");
	const int NPT=21;
	double q,q1,numout,numside,numlong,cf,N,chiout=0.0,chiside=0.0,chilong=0.0,numer,denom;
	vector<double> cfout1(NPT),cfout(NPT),errorout(NPT),cfside1(NPT),cfside(NPT),errorside(NPT),cflong1(NPT),cflong(NPT),errorlong(NPT);
	int ipt;
	numer=denom=0.0;
	for(ipt=0;ipt<NPT;ipt++){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&q,&cfout[ipt],&numout,&cfside[ipt],&numside,&cflong[ipt],&numlong, &cf, &N);
		fscanf(fptr1,"%lf %lf %lf %lf",&q1,&cfout1[ipt],&cfside1[ipt],&cflong1[ipt]);
		cfout1[ipt]-=1.0;
		cfside[ipt]-=1.0;
		cflong[ipt]-=1.0;
		errorout[ipt]=1.0/double(numout);
		errorside[ipt]=1.0/double(numside);
		errorlong[ipt]=1.0/double(numlong);
		numer+=cfout1[ipt]*cfout[ipt]/errorout[ipt]
			+cfside1[ipt]*cfside[ipt]/errorside[ipt]+cflong1[ipt]*cflong[ipt]/errorlong[ipt];
		denom+=cfout[ipt]*cfout[ipt]/errorout[ipt]
			+cfside[ipt]*cfside[ipt]/errorside[ipt]+cflong[ipt]*cflong[ipt]/errorlong[ipt];
	}
	A=numer/denom;
	for(ipt=0;ipt<NPT;ipt++){
		chiout+=pow(A*cfout[ipt]-cfout1[ipt],2)/errorout[ipt];
		chiside+=pow(A*cfside[ipt]-cfside1[ipt],2)/errorside[ipt];
		chilong+=pow(A*cflong[ipt]-cflong1[ipt],2)/errorlong[ipt];
	}
	return chiout+chiside+chilong;
}

int main(int argc,char *argv[]){
	double A, tau, R;
	if (argc != 3) {
		printf("Usage: balance_hbt tau R\n");
		exit(-1);
	}
	else{
		tau=atof(argv[1]);
		R=atof(argv[2]);
	}
	
	char scriptCommand[200];
	sprintf(scriptCommand,"runner_cf.sh %g %g",tau,R);
	system(scriptCommand);
	double chi2=CalcChiSquared_pion(A);
	printf("%g\n",chi2);
	
	return 0;	
}