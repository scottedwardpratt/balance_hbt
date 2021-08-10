#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	const int NRUNS=24,NQ=100;
	int iq,irun;
	FILE *fptr;
	char filename[120],dummy[120];
	double cout,cside,clong,denomout,denomside,denomlong;
	vector<double> q;
	vector<double> cf_out,cf_side,cf_long;
	vector<double> denom_out,denom_side,denom_long;
	q.resize(NQ);
	cf_out.resize(NQ);
	cf_side.resize(NQ);
	cf_long.resize(NQ);
	denom_out.resize(NQ);
	denom_side.resize(NQ);
	denom_long.resize(NQ);
		
	for(iq=0;iq<NQ;iq++){
		cf_out[iq]=cf_side[iq]=cf_long[iq]=denom_out[iq]=denom_side[iq]=denom_long[iq]=0.0;
	}
	for(irun=0;irun<NRUNS;irun++){
		sprintf(filename,"results/cf%d_outsidelong.dat",irun);
		printf("will read from %s\n",filename);
		fptr=fopen(filename,"r");
		fgets(dummy,120,fptr);
		for(iq=0;iq<NQ;iq++){
			fscanf(fptr,"%lf",&q[iq]);
			fscanf(fptr,"%lf %lf %lf %lf %lf %lf",
			&cout,&denomout,&cside,&denomside,&clong,&denomlong);
			if(cout!=cout)
				cout=0.0;
			if(cside!=cside)
				cside=0.0;
			if(clong!=clong)
				clong=0.0;
			cf_out[iq]+=cout*denomout;
			denom_out[iq]+=denomout;
			cf_side[iq]+=cside*denomside;
			denom_side[iq]+=denomside;
			cf_long[iq]+=clong*denomlong;
			denom_long[iq]+=denomlong;
		}
		fclose(fptr);
	}

	sprintf(filename,"results/cf_outsidelong.dat");
	fptr=fopen(filename,"w");
	for(iq=0;iq<NQ;iq++){
		cf_out[iq]=cf_out[iq]/(denom_out[iq]);
		if(cf_out[iq]!=cf_out[iq]){
			cf_out[iq]=0.0;
			printf("denom_out[%d]=%g",iq,denom_out[iq]);
		}
		cf_side[iq]=cf_side[iq]/(denom_side[iq]);
		if(cf_side[iq]!=cf_side[iq]){
			cf_side[iq]=0.0;
			printf("denom_side[%d]=%g",iq,denom_side[iq]);
		}
		cf_long[iq]=cf_long[iq]/(denom_long[iq]);
		if(cf_long[iq]!=cf_long[iq]){
			cf_long[iq]=0.0;
			printf("denom_long[%d]=%g",iq,denom_long[iq]);
		}
		fprintf(fptr,"%7.2f %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n",q[iq],
		cf_out[iq],denom_out[iq],cf_side[iq],denom_side[iq],cf_long[iq],denom_long[iq]);
	}
	fclose(fptr);
	
	return 0;
}


