#include "balhbt.h"
using namespace std;

CBF::CBF(int NYBINSset,int NPHIBINSset){
	NYBINS=NYBINSset; NPHIBINS=NPHIBINSset;
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
	YMAX=0.1*NYBINS;
	DELY=0.1;
	
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

void CBF::Increment(vector<CHBTPart *> &partvec,vector<CHBTPart *> &partvecprime,double weight){
	CHBTPart *part,*partprime;
	int ipair,iy,iphi,pid,pidprime;
	double y,yprime,phi,phiprime,dy,dphi,qqprime;
	for(ipair=0;ipair<4;ipair++){
		if(ipair==0){
			part=partvec[0]; partprime=partvecprime[0];
			pid=partvec[0]->resinfo->code; pidprime=partvecprime[0]->resinfo->code;
		}
		else if(ipair==1){
			part=partvec[0]; partprime=partvecprime[1];
			pid=partvec[0]->resinfo->code; pidprime=partvecprime[1]->resinfo->code;
		}
		else if(ipair==2){
			part=partvec[1]; partprime=partvecprime[0];
			pid=partvec[1]->resinfo->code; pidprime=partvecprime[0]->resinfo->code;
		}
		else{
			part=partvec[1]; partprime=partvecprime[1];
			pid=partvec[1]->resinfo->code; pidprime=partvecprime[1]->resinfo->code;
		}
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
	
}
