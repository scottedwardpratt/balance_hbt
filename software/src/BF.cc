#include "balhbt.h"
using namespace std;

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
	QINVMAX=DELQINV*NQINVBINS;
	Zero();
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
	picount=Kcount=pcount=Ntest=0;
	
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
	fprintf(balhbt->logfile,"picount=%lld, Kcount=%lld, pcount=%lld, Ntest=%lld\n",picount,Kcount,pcount,Ntest);
	fflush(balhbt->logfile);
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
