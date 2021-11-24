#include "balhbt.h"
//using namespace std;

CCF_Arrays::CCF_Arrays(int NYset,double DELYset,int NPHIset,int NQset,double DELQset){
	NY=NYset;
	DELY=DELYset;
	NPHI=NPHIset;
	NQ=NQset;
	DELQ=DELQset;
	DELPHI=180.0/double(NPHI);
	cf_y.resize(NY,0.0); denom_y.resize(NY,0.0);
	cf_phi.resize(NPHI,0.0); denom_phi.resize(NPHI,0.0);
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
	}
	iphi=lrint(floor(delphi/DELPHI));
	if(iphi<NPHI){
		cf_phi[iphi]+=weight;
		denom_phi[iphi]+=1.0;
	}
}

void CCF_Arrays::Print(){
	printf("-----------\n");
	printf("DELQ=%g, NQ=%d, DELY=%g, NY=%d, NPHI=%d\n",DELQ,NQ,DELY,NY,NPHI);
	printf("-----------\n");
}

void CCF_Arrays::WriteResultsCF(string dirname,long long int denom_count,int run_number){
	FILE *fptr;
	double delq,dely,delphi;
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

		fprintf(fptr,"%7.4f %12.9f %12.0f\n",dely,cf_y[iy]/denom_y[iy],denom_y[iy]);
	}
	fclose(fptr);
	filename=dirname+"cf"+to_string(run_number)+"_phi.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   delphi     C(delphi)    Nphi\n");
	for(int iphi=0;iphi<NPHI;iphi++){
		delphi=(iphi+0.5)*DELPHI;
		fprintf(fptr,"%7.4f %12.9f %12.0f\n",delphi,cf_phi[iphi]/denom_phi[iphi],denom_phi[iphi]);
	}
	fclose(fptr);
}

void CCF_Arrays::WriteResultsBF(string dirname,long long int denom_count,int run_number){
	FILE *fptr;
	double delq,dely,delphi,denom,norm=0.0;
	string filename;
	string command;
	command="mkdir -p "+dirname;
	system(command.c_str());
	filename=dirname+"bf"+to_string(run_number)+"_outsidelong.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   q     Cout       Nout      Cside       Nside       Clong        Nlong\n");
	for(int iq=0;iq<NQ;iq++){
		delq=(0.5+iq)*DELQ;
		denom=0.5*denom_count*DELPHI;
		fprintf(fptr,"%7.3f %12.9f %12.0f %12.9f %12.0f %12.9f %12.0f %12.9f %12.0f\n",
		delq,
		cf_out[iq]/denom,denom_out[iq],cf_side[iq]/denom,denom_side[iq],
		cf_long[iq]/denom,denom_long[iq],cf_inv[iq]/denom,denom_inv[iq]);
	}
	fclose(fptr);
	filename=dirname+"bf"+to_string(run_number)+"_y.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   dely     C(y)    Ny\n");
	for(int iy=0;iy<NY;iy++){
		dely=(iy+0.5)*DELY;
		denom=0.5*denom_count*DELY;
		fprintf(fptr,"%7.4f %12.9f %12.0f\n",dely,cf_y[iy]/denom,denom_y[iy]);
		norm+=DELY*cf_y[iy]/denom;
	}
	printf("%s: norm=%g\n",dirname.c_str(),norm);
	fclose(fptr);
	filename=dirname+"bf"+to_string(run_number)+"_phi.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#   delphi     C(delphi)    Nphi\n");
	for(int iphi=0;iphi<NPHI;iphi++){
		delphi=(iphi+0.5)*DELPHI;
		denom=0.5*denom_count*DELPHI*PI/180.0;
		fprintf(fptr,"%7.4f %12.9f %12.0f\n",delphi,cf_phi[iphi]/denom,denom_phi[iphi]);
	}
	fclose(fptr);
}
