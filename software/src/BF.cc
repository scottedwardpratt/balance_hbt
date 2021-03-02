#include "balhbt.h"
using namespace std;

CBF::CBF(CparameterMap &parmap){
	NYBINS=parmap.getD("BF_NYBINS",20);
	NPHIBINS=parmap.getD("BF_NPHIBINS",18);
	DELY=parmap.getD("BF_DELY",0.1);
	BFy_pipi.resize(NYBINS);
	BFy_piK.resize(NYBINS);
	BFy_pip.resize(NYBINS);
	BFy_KK.resize(NYBINS);
	BFy_Kp.resize(NYBINS);
	BFy_pp.resize(NYBINS);
	
	BFphi_pipi.resize(NYBINS);
	BFphi_piK.resize(NYBINS);
	BFphi_pip.resize(NYBINS);
	BFphi_KK.resize(NYBINS);
	BFphi_Kp.resize(NYBINS);
	BFphi_pp.resize(NYBINS);

	DENOMy_pipi.resize(NYBINS);
	DENOMy_piK.resize(NYBINS);
	DENOMy_pip.resize(NYBINS);
	DENOMy_KK.resize(NYBINS);
	DENOMy_Kp.resize(NYBINS);
	DENOMy_pp.resize(NYBINS);
	
	DENOMphi_pipi.resize(NYBINS);
	DENOMphi_piK.resize(NYBINS);
	DENOMphi_pip.resize(NYBINS);
	DENOMphi_KK.resize(NYBINS);
	DENOMphi_Kp.resize(NYBINS);
	DENOMphi_pp.resize(NYBINS);
	
	DELPHI=PI/double(NPHIBINS);
	YMAX=DELY*NYBINS;
	
	Zero();
}

void CBF::Zero(){
	
	BFy_pipi.assign(NYBINS,0.0);
	BFy_piK.assign(NYBINS,0.0);
	BFy_pip.assign(NYBINS,0.0);
	BFy_KK.assign(NYBINS,0.0);
	BFy_Kp.assign(NYBINS,0.0);
	BFy_pp.assign(NYBINS,0.0);
	
	BFphi_pipi.assign(NYBINS,0.0);
	BFphi_piK.assign(NYBINS,0.0);
	BFphi_pip.assign(NYBINS,0.0);
	BFphi_KK.assign(NYBINS,0.0);
	BFphi_Kp.assign(NYBINS,0.0);
	BFphi_pp.assign(NYBINS,0.0);
	
	DENOMy_pipi.assign(NYBINS,0.0);
	DENOMy_piK.assign(NYBINS,0.0);
	DENOMy_pip.assign(NYBINS,0.0);
	DENOMy_KK.assign(NYBINS,0.0);
	DENOMy_Kp.assign(NYBINS,0.0);
	DENOMy_pp.assign(NYBINS,0.0);
	
	DENOMphi_pipi.assign(NYBINS,0.0);
	DENOMphi_piK.assign(NYBINS,0.0);
	DENOMphi_pip.assign(NYBINS,0.0);
	DENOMphi_KK.assign(NYBINS,0.0);
	DENOMphi_Kp.assign(NYBINS,0.0);
	DENOMphi_pp.assign(NYBINS,0.0);
	
}

void CBF::Evaluate(vector<CHBTPart *> &partvec,vector<CHBTPart *> &productvec,
vector<CHBTPart *> &partprimevec,vector<CHBTPart *> &productprimevec,double balweight){
	double dely,dely0,delytarget,hbtweight;
	unsigned int i,iprime,j,k;
	CHBTPart *part,*partprime;
	for(i=0;i<productvec.size();i++){
		part=productvec[i];
		for(iprime=0;iprime<productprimevec.size();iprime++){
			partprime=productprimevec[iprime];
			dely0=part.GetRapidity()-partprime.GetRapidity();
			delytarget=YMAX*(1.0-2.0*randy->ran());
			dely=delytarget-dely0;
			for(j=0;j<2;j++){
				partvec[j].BjBoost(dely);
			}
			for(k=0;k<productvec.size();k++){
				productvec[k].BjBoost(dely);
			}
			hbtweight=hbtcalc->GetPsiSquared(partvec[0],partprimevec[0])
				*hbtcalc->GetPsiSquared(partvec[0],partprimevec[1])
					hbtcalc->GetPsiSquared(partvec[1],partprimevec[0])
						hbtcalc->GetPsiSquared(partvec[1],partprimevec[1]);
			Increment(part,partprime,balweight*hbtweight);
		}
	}
}

void CBF::Increment(CHBTPart *part,CHBTPart partprime,double weight){
	int iy,iphi,pid,pidprime;
	double y,yprime,phi,phiprime,dy,dphi,qqprime;
	pid=part->resinfo->code;
	pidprime=partprime->resinfo->code;
	pid=abs(pid);
	pidprime=abs(pidprime);
	if((pid==211 || pid==321 || pid==2212) && (pidprime==211 || pidprime==321 || pidprime==2212)){
		qqprime=part->resinfo->charge*partprime->resinfo->charge;
		y=atanh(part->p[3]/part->p[0]);
		phi=atan2(part->p[2],part->p[1]);
		yprime=atanh(partprime->p[3]/partprime->p[0]);
		phiprime=atan2(partprime->p[2],partprime->p[1]);
		dy=y-yprime;
		if(dy<-YMAX)
			dy+=2.0*YMAX;
		if(dy>YMAX)
			dy-=2.0*YMAX;
		dy=fabs(dy);
		iy=lrint(floor(dy/DELY));
		if(iy<NYBINS){
			dphi=phi-phiprime;
			if(dphi<-PI)
				dphi+=2.0*PI;
			if(dphi>PI)
				dphi-=2.0*PI;
			dphi=fabs(dphi);
			iphi=lrint(floor((dphi*180.0/PI)/DELPHI));
			if(pid==211 && pidprime==211){
				BFy_pipi[iy]-=qqprime*weight;
				BFphi_pipi[iphi]-=qqprime*weight;
				DENOMy_pipi[iy]+=1.0;
				DENOMphi_pipi[iphi]+=1.0;
			}
			else if((pid==211 && pidprime==321) || (pid==321 && pidprime==211)){
				BFy_piK[iy]-=qqprime*weight;
				BFphi_piK[iphi]-=qqprime*weight;
				DENOMy_piK[iy]+=1.0;
				DENOMphi_piK[iphi]+=1.0;
			}
			else if((pid==211 && pidprime==2212) || (pid==2212 && pidprime==211)){
				BFy_pip[iy]-=qqprime*weight;
				BFphi_pip[iphi]-=qqprime*weight;
				DENOMy_pipi[iy]+=1.0;
				DENOMphi_pip[iphi]+=1.0;
			}
			else if(pid==321 && pidprime==321){
				BFy_KK[iy]-=qqprime*weight;
				BFphi_KK[iphi]-=qqprime*weight;
				DENOMy_KK[iy]+=1.0;
				DENOMphi_KK[iphi]+=1.0;
			}
			else if((pid==321 && pidprime==2212) || (pid==2212 && pidprime==321)){
				BFy_Kp[iy]-=qqprime*weight;
				BFphi_Kp[iphi]-=qqprime*weight;
				DENOMy_Kp[iy]+=1.0;
				DENOMphi_Kp[iphi]+=1.0;
			}
			else if(pid==2212 && pidprime==2212){
				BFy_pp[iy]-=qqprime*weight;
				BFphi_pp[iphi]-=qqprime*weight;
				DENOMy_pp[iy]+=1.0;
				DENOMphi_pp[iphi]+=1.0;
			}
		}
	}	
}
