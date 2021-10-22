#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
//#include <conio.h>
#include "parametermap.h"
#include <stdlib.h>


double  CalcChiSquared_pion(double A){
	FILE *fptr=fopen("results/pipluspiplus/cf_outsidelong.dat","r");
	FILE *fptr1=fopen("cf_exp.dat","r");
	double q,q1,cfout,cfout1,cfside1,cflong1,numout,cfside,numside,cflong,numlong,errorout, errorside, errorlong, cf, N, chiout=0.0,chiside=0.0,chilong=0.0;
	int ipt;
	for(ipt=1;ipt<22;ipt++){
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&q,&cfout,&numout,&cfside,&numside,&cflong,&numlong, &cf, &N);
		fscanf(fptr1,"%lf %lf %lf %lf",&q1,&cfout1,&cfside1,&cflong1);
		errorout=1.0/sqrt(numout);
		errorside=1.0/sqrt(numside);
		errorlong=1.0/sqrt(numlong);
		
		chiout+=pow(1+A*cfout-cfout1,2)/errorout;
		chiside+=pow(1+A*cfside-cfside1,2)/errorside;
		chilong+=pow(1+A*cflong-cflong1,2)/errorlong;
		
	}
	return chiout+chiside+chilong;
}

int main(int argc,char *argv[]){
	double A, tau, R;
	if (argc != 4) {
		printf("Usage: balance_hbt A tau R\n");
		exit(-1);
	}
	else{
		A=atof(argv[1]);
		tau=atof(argv[2]);
		R=atof(argv[3]);
	}
//	CparameterMap parmap;
//	parmap.ReadParsFromFile("parameters/bwpars.txt");
//	parmap.set("BW_TAU",tau);
//	parmap.set("BW_RPERP",R);
	
    string scriptCommand = "runner_cf.sh";

    // Loop through and add all the arguments.
    for (int i = 2; i < argc; i++){
        scriptCommand = scriptCommand + " " + argv[i];	
    }	
	system(scriptCommand.c_str());
	double chi2=CalcChiSquared_pion(A);
	printf("%04.2f\n",chi2);
	
	return 0;	
}