#include "balhbt.h"
using namespace std;

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
		cf_y[iq]=denom_y[iq]=0.0;
		cf_phi[iq]=denom_phi[iq]=0.0;
		cf_inv[iq]=denom_inv[iq]=0.0;
		cf_out[iq]=denom_out[iq]=0.0;
		cf_side[iq]=denom_side[iq]=0.0;
		cf_long[iq]=denom_long[iq]=0.0;
	}
}

void CCF_Arrays::Increment(double dely,double delphi,double qinv,double qout,double qside,double qlong,double weight){
	double QDIRCUT=15.0,qother;
	int iq,iy,iphi;
	iq=lrint(floor(qinv/DELQ));
	if(iq<NQ){
		cf_inv[iq]+=weight;
		denom_inv[iq]+=1.0;
	}
	qother=sqrt(qside*qside+qlong*qlong);
	if(qother<QDIRCUT){
		iq=lrint(floor(qout/DELQ));
		if(iq<NQ){
			cf_out[iq]+=weight;
			denom_out[iq]+=1.0;
		}
	}
	qother=sqrt(qout*qout+qlong*qlong);
	if(qother<QDIRCUT){
		iq=lrint(floor(qside/DELQ));
		if(iq<NQ){
			cf_side[iq]+=weight;
			denom_side[iq]+=1.0;
		}
	}
	qother=sqrt(qside*qside+qout*qout);
	if(qother<QDIRCUT){
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

CBF::CBF(CBalHBT *balhbtset){
	balhbt=balhbtset;
	CBF_init(balhbt->parmap);
	CHEAPPSISQUARED=parmap->getB("BF_CHEAPPSISQUARED",false);
	UseAllWFsForCF=parmap->getB("BF_USEALLWFSFORCF",true);
}

void CBF::CBF_init(CparameterMap *parmapin){
	parmap=parmapin;

	NYBINS=parmap->getI("BF_NYBINS",20);
	NPHIBINS=parmap->getI("BF_NPHIBINS",18);
	DELY=parmap->getD("BF_DELY",0.1);
	DELQINV=parmap->getD("BF_DELQINV",5.0);
 	NQINVBINS=parmap->getI("BF_NQINVBINS",100);
	acceptancebal=new CAcceptanceBal(parmap);
	acceptancebal->balhbt=balhbt;
	
	CF_pipluspiplus=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_pipluspiminus=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_piplusKplus=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_piplusKminus=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_piplusp=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_pipluspbar=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_KplusKplus=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_KplusKminus=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_Kplusp=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_Kpluspbar=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_pp=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	CF_ppbar=new CCF_Arrays(NYBINS,DELY,NPHIBINS,NQINVBINS,DELQINV);
	
	DELPHI=180.0/double(NPHIBINS);
	YMAX=DELY*NYBINS;
	QINVMAX=DELQINV*NQINVBINS;

}

void CBF::Zero(){
	
	CF_pipluspiplus->Zero();
	CF_pipluspiminus->Zero();
	CF_piplusKplus->Zero();
	CF_piplusKminus->Zero();
	CF_piplusp->Zero();
	CF_pipluspbar->Zero();
	CF_KplusKplus->Zero();
	CF_KplusKminus->Zero();
	CF_Kplusp->Zero();
	CF_Kpluspbar->Zero();
	CF_pp->Zero();
	CF_ppbar->Zero();
	
	picount=Kcount=pcount=0;
	
}

void CBF::WriteResults(int run_number){
	string dirname;
	string type;
	type="pipluspiplus"; dirname="results/"+type+"/";
	CF_pipluspiplus->WriteResults(dirname,run_number);
	type="pipluspiminus"; dirname="results/"+type+"/";
	CF_pipluspiminus->WriteResults(dirname,run_number);
	type="piplusKplus"; dirname="results/"+type+"/";
	CF_piplusKplus->WriteResults(dirname,run_number);
	type="piplusKminus"; dirname="results/"+type+"/";
	CF_piplusKminus->WriteResults(dirname,run_number);
	type="piplusp"; dirname="results/"+type+"/";
	CF_piplusp->WriteResults(dirname,run_number);
	type="pipluspbar"; dirname="results/"+type+"/";
	CF_pipluspbar->WriteResults(dirname,run_number);
	type="KplusKplus"; dirname="results/"+type+"/";
	CF_KplusKplus->WriteResults(dirname,run_number);
	type="KplusKminus"; dirname="results/"+type+"/";
	CF_KplusKminus->WriteResults(dirname,run_number);
	type="Kplusp"; dirname="results/"+type+"/";
	CF_Kplusp->WriteResults(dirname,run_number);
	type="Kpluspbar"; dirname="results/"+type+"/";
	CF_Kpluspbar->WriteResults(dirname,run_number);
	type="pp"; dirname="results/"+type+"/";
	CF_pp->WriteResults(dirname,run_number);
	type="ppbar"; dirname="results/"+type+"/";
	CF_ppbar->WriteResults(dirname,run_number);
}

void CCF_Arrays::WriteResults(string dirname,int run_number){
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

double CBF::Getqinv(CHBTPart *part,CHBTPart *partprime){
	double P2,Pdotq,Qinv2,m=part->resinfo->mass,mprime=partprime->resinfo->mass;
	int alpha;
	FourVector P,q;
	for(alpha=0;alpha<4;alpha++){
		P[alpha]=part->p[alpha]+partprime->p[alpha];
		q[alpha]=part->p[alpha]-partprime->p[alpha];
	}
	P2=P[0]*P[0]-P[1]*P[1]-P[2]*P[2]-P[3]*P[3];
	Pdotq=m*m-mprime*mprime;
	Qinv2=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
	Qinv2=Qinv2-Pdotq*Pdotq/P2;
	return 0.5*sqrt(-Qinv2);
}

double CBF::CheapPsiSquared(CHBTPart *part,CHBTPart *partprime){
	double answer,qinv,Rinv=20.0;
	if(part->resinfo->code==partprime->resinfo->code && abs(part->resinfo->code)==211){
		qinv=Getqinv(part,partprime);
		if(qinv<2000){
			answer=1.0+exp(-qinv*qinv*Rinv*Rinv/(HBARC*HBARC));
		}
		else answer=1.0;
	}
	else
		answer=1.0;
	return answer;
}
