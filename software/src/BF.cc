#include "balhbt.h"
using namespace std;

double CBF::netweight=0.0;

CBF::CBF(CparameterMap *parmapin){
	parmap=parmapin;

	NYBINS=parmap->getI("BF_NYBINS",20);
	NPHIBINS=parmap->getI("BF_NPHIBINS",18);
	DELY=parmap->getD("BF_DELY",0.1);
	DELQINV=parmap->getD("BF_DELQINV",5.0);
	NQINVBINS=parmap->getI("BF_NQINVBINS",100);
	acceptancebal=new CAcceptanceBal(parmap);
	
	BFy_pipi.resize(NYBINS);
	BFy_piK.resize(NYBINS);
	BFy_pip.resize(NYBINS);
	BFy_KK.resize(NYBINS);
	BFy_Kp.resize(NYBINS);
	BFy_pp.resize(NYBINS);
	
	BFphi_pipi.resize(NPHIBINS);
	BFphi_piK.resize(NPHIBINS);
	BFphi_pip.resize(NPHIBINS);
	BFphi_KK.resize(NPHIBINS);
	BFphi_Kp.resize(NPHIBINS);
	BFphi_pp.resize(NPHIBINS);
	
	BFqinv_pipi.resize(NQINVBINS);
	BFqinv_piK.resize(NQINVBINS);
	BFqinv_pip.resize(NQINVBINS);
	BFqinv_KK.resize(NQINVBINS);
	BFqinv_Kp.resize(NQINVBINS);
	BFqinv_pp.resize(NQINVBINS);

	DENOMy_pipi.resize(NYBINS);
	DENOMy_piK.resize(NYBINS);
	DENOMy_pip.resize(NYBINS);
	DENOMy_KK.resize(NYBINS);
	DENOMy_Kp.resize(NYBINS);
	DENOMy_pp.resize(NYBINS);
	
	DENOMphi_pipi.resize(NPHIBINS);
	DENOMphi_piK.resize(NPHIBINS);
	DENOMphi_pip.resize(NPHIBINS);
	DENOMphi_KK.resize(NPHIBINS);
	DENOMphi_Kp.resize(NPHIBINS);
	DENOMphi_pp.resize(NPHIBINS);
	
	DENOMqinv_pipi.resize(NQINVBINS);
	DENOMqinv_piK.resize(NQINVBINS);
	DENOMqinv_pip.resize(NQINVBINS);
	DENOMqinv_KK.resize(NQINVBINS);
	DENOMqinv_Kp.resize(NQINVBINS);
	DENOMqinv_pp.resize(NQINVBINS);
	
	CFqinv_pipluspiplus.resize(NQINVBINS);
	CFqinv_pipluspiminus.resize(NQINVBINS);
	CFqinv_KplusKplus.resize(NQINVBINS);
	CFqinv_KplusKminus.resize(NQINVBINS);
	CFqinv_pp.resize(NQINVBINS);
	CFqinv_ppbar.resize(NQINVBINS);
	CF_DENOMqinv_pipluspiplus.resize(NQINVBINS);
	CF_DENOMqinv_pipluspiminus.resize(NQINVBINS);
	CF_DENOMqinv_KplusKplus.resize(NQINVBINS);
	CF_DENOMqinv_KplusKminus.resize(NQINVBINS);
	CF_DENOMqinv_pp.resize(NQINVBINS);
	CF_DENOMqinv_ppbar.resize(NQINVBINS);
	
	DELPHI=PI/double(NPHIBINS);
	YMAX=DELY*NYBINS;
	QINVMAX=DELQINV*NQINVBINS;
	
	Zero();
}

void CBF::Zero(){
	
	BFy_pipi.assign(NYBINS,0.0);
	BFy_piK.assign(NYBINS,0.0);
	BFy_pip.assign(NYBINS,0.0);
	BFy_KK.assign(NYBINS,0.0);
	BFy_Kp.assign(NYBINS,0.0);
	BFy_pp.assign(NYBINS,0.0);
	
	BFphi_pipi.assign(NPHIBINS,0.0);
	BFphi_piK.assign(NPHIBINS,0.0);
	BFphi_pip.assign(NPHIBINS,0.0);
	BFphi_KK.assign(NPHIBINS,0.0);
	BFphi_Kp.assign(NPHIBINS,0.0);
	BFphi_pp.assign(NPHIBINS,0.0);
	
	DENOMy_pipi.assign(NYBINS,0.0);
	DENOMy_piK.assign(NYBINS,0.0);
	DENOMy_pip.assign(NYBINS,0.0);
	DENOMy_KK.assign(NYBINS,0.0);
	DENOMy_Kp.assign(NYBINS,0.0);
	DENOMy_pp.assign(NYBINS,0.0);
	
	DENOMphi_pipi.assign(NPHIBINS,0.0);
	DENOMphi_piK.assign(NPHIBINS,0.0);
	DENOMphi_pip.assign(NPHIBINS,0.0);
	DENOMphi_KK.assign(NPHIBINS,0.0);
	DENOMphi_Kp.assign(NPHIBINS,0.0);
	DENOMphi_pp.assign(NPHIBINS,0.0);
	
	picount=Kcount=pcount=0;
	
}

void CBF::WriteResults(int run_number){
	char filename[100];
	sprintf(filename,"results/bf%d_y.dat",run_number);
	FILE *fptr=fopen(filename,"w");
	for(int iy=0;iy<NYBINS;iy++){
		fprintf(fptr,"%7.3f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(iy+0.5)*DELY,
		BFy_pipi[iy],DENOMy_pipi[iy],BFy_piK[iy],DENOMy_piK[iy],BFy_pip[iy],DENOMy_pip[iy],
		BFy_KK[iy],DENOMy_KK[iy],BFy_Kp[iy],DENOMy_Kp[iy],BFy_pp[iy],DENOMy_pp[iy]);
	}
	fclose(fptr);
	sprintf(filename,"results/bf%d_phi.dat",run_number);
	fptr=fopen(filename,"w");
	for(int iphi=0;iphi<NPHIBINS;iphi++){
		fprintf(fptr,"%7.3f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(iphi+0.5)*DELPHI*180.0/PI,
		BFphi_pipi[iphi],DENOMphi_pipi[iphi],BFphi_piK[iphi],DENOMphi_piK[iphi],BFphi_pip[iphi],DENOMphi_pip[iphi],
		BFphi_KK[iphi],DENOMphi_KK[iphi],BFphi_Kp[iphi],DENOMphi_Kp[iphi],BFphi_pp[iphi],DENOMphi_pp[iphi]);
	}
	fclose(fptr);
	sprintf(filename,"results/bf%d_qinv.dat",run_number);
	fptr=fopen(filename,"w");
	for(int iqinv=0;iqinv<NQINVBINS;iqinv++){
		fprintf(fptr,"%7.3f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",(iqinv+0.5)*DELQINV,
		BFqinv_pipi[iqinv],DENOMqinv_pipi[iqinv],BFqinv_piK[iqinv],DENOMqinv_piK[iqinv],BFqinv_pip[iqinv],DENOMqinv_pip[iqinv],
		BFqinv_KK[iqinv],DENOMqinv_KK[iqinv],BFqinv_Kp[iqinv],DENOMqinv_Kp[iqinv],BFqinv_pp[iqinv],DENOMqinv_pp[iqinv]);
	}
	fclose(fptr);
	sprintf(filename,"results/bf%d_hadcount.dat",run_number);
	fptr=fopen(filename,"w");
	fprintf(fptr,"%lld  %lld  %lld\n",picount,Kcount,pcount);
	fclose(fptr);
	
	sprintf(filename,"results/cf%d_qinv.dat",run_number);
	fptr=fopen(filename,"w");
	for(int iqinv=0;iqinv<NQINVBINS;iqinv++){
		fprintf(fptr,"%7.3f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",
		(iqinv+0.5)*DELQINV,
		CFqinv_pipluspiplus[iqinv]/CF_DENOMqinv_pipluspiplus[iqinv],CF_DENOMqinv_pipluspiplus[iqinv],
		CFqinv_pipluspiminus[iqinv]/CF_DENOMqinv_pipluspiminus[iqinv],CF_DENOMqinv_pipluspiminus[iqinv],
		CFqinv_KplusKplus[iqinv]/CF_DENOMqinv_KplusKplus[iqinv],CF_DENOMqinv_KplusKplus[iqinv],
		CFqinv_KplusKminus[iqinv]/CF_DENOMqinv_KplusKminus[iqinv],CF_DENOMqinv_KplusKminus[iqinv],
		CFqinv_pp[iqinv]/CF_DENOMqinv_pp[iqinv],CF_DENOMqinv_pp[iqinv],
		CFqinv_ppbar[iqinv]/CF_DENOMqinv_ppbar[iqinv],CF_DENOMqinv_ppbar[iqinv]);		
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
	double answer,qinv,Rinv=4.0;
	if(part->resinfo->code==partprime->resinfo->code){
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
