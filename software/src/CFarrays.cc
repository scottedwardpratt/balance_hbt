#include "balhbt.h"
//using namespace std;

CCF_Arrays::CCF_Arrays(int NYset,double DELYset,int NPHIset,int NQset,double DELQset){
	NY=NYset;
	DELY=DELYset;
	NPHI=NPHIset;
	NQ=NQset;
	DELQ=DELQset;
	DELPHI=180.0/double(NPHI);
	cf_y.resize(NY,0.0); denom_y.resize(NY,0.0); error_y.resize(NY,0.0);
	cf_phi.resize(NPHI,0.0); denom_phi.resize(NPHI,0.0); error_phi.resize(NPHI,0.0);
	cf_inv.resize(NQ,0.0); denom_inv.resize(NQ,0.0);
	cf_out.resize(NQ,0.0); denom_out.resize(NQ,0.0);
	cf_side.resize(NQ,0.0); denom_side.resize(NQ,0.0);
	cf_long.resize(NQ,0.0); denom_long.resize(NQ,0.0);
}

void CCF_Arrays::Zero(){
	for(int iq=0;iq<NQ;iq++){
		cf_inv[iq]=denom_inv[iq]=0.0;
		cf_out[iq]=denom_out[iq]=0.0;
		cf_side[iq]=denom_side[iq]=0.0;
		cf_long[iq]=denom_long[iq]=0.0;
	}
	for(int iy=0;iy<NY;iy++){
		cf_y[iy]=denom_y[iy]=0.0;
	}
	for(int iphi=0;iphi<NPHI;iphi++){
		cf_phi[iphi]=denom_phi[iphi]=0.0;
	}
}

void CCF_Arrays::Increment(double dely,double delphi,double qinv,double qout,double qside,double qlong,double weight){
	double QDIRCUT=15.0;
	int iq,iy,iphi;
	iq=lrint(floor(qinv/DELQ));
	if(iq<NQ){
		cf_inv[iq]+=weight;
		denom_inv[iq]+=1.0;
	}
	if(fabs(qside)<QDIRCUT && fabs(qlong)<QDIRCUT){
		iq=lrint(floor(qout/DELQ));
		if(iq<NQ){
			cf_out[iq]+=weight;
			denom_out[iq]+=1.0;
		}
	}
	if(fabs(qout)<QDIRCUT && fabs(qlong)<QDIRCUT){
		iq=lrint(floor(qside/DELQ));
		if(iq<NQ){
			cf_side[iq]+=weight;
			denom_side[iq]+=1.0;
		}
	}
	if(fabs(qout)<QDIRCUT && fabs(qside)<QDIRCUT){
		iq=lrint(floor(qlong/DELQ));
		if(iq<NQ){
			cf_long[iq]+=weight;
			denom_long[iq]+=1.0;
		}
	}
	iy=lrint(floor(dely/DELY));
	if(iy<NY){
		cf_y[iy]+=weight;
		denom_y[iy]+=1.0;
		error_y[iy]+=weight*weight;
	}
	iphi=lrint(floor(delphi/DELPHI));
	if(iphi<NPHI){
		cf_phi[iphi]+=weight;
		denom_phi[iphi]+=1.0;
		error_phi[iphi]+=weight*weight;
	}
}

void CCF_Arrays::Print(){
	printf("-----------\n");
	printf("DELQ=%g, NQ=%d, DELY=%g, NY=%d, NPHI=%d\n",DELQ,NQ,DELY,NY,NPHI);
	printf("-----------\n");
}

void CCF_Arrays::WriteResults(string dirname,int run_number){
	FILE *fptr;
	double delq,dely,delphi,sigma;
	string filename;
	string command;
	command="mkdir -p "+dirname;
	system(command.c_str());
	filename=dirname+"cf"+to_string(run_number)+"_outsidelong.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   q     Cout       Nout      Cside       Nside       Clong        Nlong\n");
	for(int iq=0;iq<NQ;iq++){
		delq=(0.5+iq)*DELQ;
		fprintf(fptr,"%7.3f %12.9f %12.0f %12.9f %12.0f %12.9f %12.0f %12.9f %12.0f\n",
		delq,
		cf_out[iq]/denom_out[iq],denom_out[iq],cf_side[iq]/denom_side[iq],denom_side[iq],
		cf_long[iq]/denom_long[iq],denom_long[iq],cf_inv[iq]/denom_inv[iq],denom_inv[iq]);
	}
	fclose(fptr);
	filename=dirname+"cf"+to_string(run_number)+"_y.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   dely     C(y)    Ny\n");
	for(int iy=0;iy<NY;iy++){
		dely=(iy+0.5)*DELY;
		sigma=error_y[iy]/denom_y[iy]-pow(cf_y[iy]/denom_y[iy],2);
		sigma=sqrt(sigma/denom_y[iy]);
		fprintf(fptr,"%7.4f %12.9f %12.0f %12.9f\n",dely,cf_y[iy]/denom_y[iy],denom_y[iy],sigma);
	}
	fclose(fptr);
	filename=dirname+"cf"+to_string(run_number)+"_phi.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   delphi     C(delphi)    Nphi\n");
	for(int iphi=0;iphi<NPHI;iphi++){
		delphi=(iphi+0.5)*DELPHI;
		sigma=error_phi[iphi]/denom_phi[iphi]-pow(cf_phi[iphi]/denom_phi[iphi],2);
		sigma=sqrt(sigma/denom_phi[iphi]);
		fprintf(fptr,"%7.4f %12.9f %12.0f %12.9f\n",delphi,cf_phi[iphi]/denom_phi[iphi],denom_phi[iphi],sigma);
	}
	fclose(fptr);
}
