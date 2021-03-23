#include "balhbt.h"
using namespace std;

CBF::CBF(CparameterMap *parmapin){
	parmap=parmapin;

	NYBINS=parmap->getI("BF_NYBINS",20);
	NPHIBINS=parmap->getI("BF_NPHIBINS",18);
	DELY=parmap->getD("BF_DELY",0.1);
	DELQINV=parmap->getD("BF_DELQINV",5.0);
	NQINVBINS=parmap->getI("BF_NQINVBINS",100);
	
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
	
}

void CBF::PartAntipart(CHBTPart *part){
	CResInfo *resinfo=part->resinfo;
	if(resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0){
		part->resinfo=part->resinfo->reslist->GetResInfoPtr(-resinfo->code);
	}
}

void CBF::Evaluate(vector<CHBTPart *> &partvec,vector<vector<CHBTPart *>> &productvec,
vector<CHBTPart *> &partprimevec,vector<vector<CHBTPart *>> &productprimevec,double balweight){
	double hbtweight,psisquared1,psisquared2,psisquared3,psisquared4;
	unsigned int i,iprime,iprod,iprodprime,jpartantipart;
	CHBTPart *part,*partprime;

	for(jpartantipart=0;jpartantipart<4;jpartantipart++){
		if(jpartantipart==2){
			for(i=0;i<2;i++){
				PartAntipart(partvec[i]);
				for(iprod=0;iprod<productvec[i].size();iprod++){
					PartAntipart(productvec[i][iprod]);
				}
			}
		}
		if(jpartantipart==1 || jpartantipart==3){
			for(iprime=0;iprime<2;iprime++){
				PartAntipart(partprimevec[iprime]);
				for(iprodprime=0;iprodprime<productprimevec[iprime].size();iprodprime++){
					PartAntipart(productprimevec[iprime][iprodprime]);
				}
			}
		}
		/*
		psisquared1=hbtcalc->GetPsiSquared(partvec[0],partprimevec[0]);
		psisquared2=hbtcalc->GetPsiSquared(partvec[0],partprimevec[1]);
		psisquared3=hbtcalc->GetPsiSquared(partvec[1],partprimevec[0]);
		psisquared4=hbtcalc->GetPsiSquared(partvec[1],partprimevec[1]);
		*/
		psisquared1=CheapPsiSquared(partvec[0],partprimevec[0]);
		psisquared2=CheapPsiSquared(partvec[0],partprimevec[1]);
		psisquared3=CheapPsiSquared(partvec[1],partprimevec[0]);
		psisquared4=CheapPsiSquared(partvec[1],partprimevec[1]);
		
		hbtweight=psisquared1*psisquared2*psisquared3*psisquared4;
		
		for(i=0;i<2;i++){
			for(iprime=0;iprime<2;iprime++){
				for(iprod=0;iprod<productvec[i].size();iprod++){
					part=productvec[i][iprod];
					for(iprodprime=0;iprodprime<productprimevec[iprime].size();iprodprime++){
						partprime=productprimevec[iprime][iprodprime];
						Increment(part,partprime,balweight*hbtweight);
					}
				}
			}
		}
	}
}

void CBF::Increment(CHBTPart *part,CHBTPart *partprime,double weight){
	int iy,iphi,iqinv,pid,pidprime;
	double y,yprime,phi,phiprime,dy,dphi,qqprime,qinv;
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
			iphi=lrint(floor(dphi/DELPHI));
			qinv=Getqinv(part,partprime);
			iqinv=lrint(floor(qinv)/DELQINV);
			if(pid==211 && pidprime==211){
				BFy_pipi[iy]-=qqprime*weight;
				BFphi_pipi[iphi]-=qqprime*weight;
				DENOMy_pipi[iy]+=1.0;
				DENOMphi_pipi[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_pipi[iqinv]-=qqprime*weight;
					DENOMqinv_pipi[iqinv]+=1.0;
				}
			}
			else if((pid==211 && pidprime==321) || (pid==321 && pidprime==211)){
				BFy_piK[iy]-=qqprime*weight;
				BFphi_piK[iphi]-=qqprime*weight;
				DENOMy_piK[iy]+=1.0;
				DENOMphi_piK[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_piK[iqinv]-=qqprime*weight;
					DENOMqinv_piK[iqinv]+=1.0;
				}
			}
			else if((pid==211 && pidprime==2212) || (pid==2212 && pidprime==211)){
				BFy_pip[iy]-=qqprime*weight;
				BFphi_pip[iphi]-=qqprime*weight;
				DENOMy_pip[iy]+=1.0;
				DENOMphi_pip[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_pip[iqinv]-=qqprime*weight;
					DENOMqinv_pip[iqinv]+=1.0;
				}
			}
			else if(pid==321 && pidprime==321){
				BFy_KK[iy]-=qqprime*weight;
				BFphi_KK[iphi]-=qqprime*weight;
				DENOMy_KK[iy]+=1.0;
				DENOMphi_KK[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_KK[iqinv]-=qqprime*weight;
					DENOMqinv_KK[iqinv]+=1.0;
				}
			}
			else if((pid==321 && pidprime==2212) || (pid==2212 && pidprime==321)){
				BFy_Kp[iy]-=qqprime*weight;
				BFphi_Kp[iphi]-=qqprime*weight;
				DENOMy_Kp[iy]+=1.0;
				DENOMphi_Kp[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_Kp[iqinv]-=qqprime*weight;
					DENOMqinv_Kp[iqinv]+=1.0;
				}
			}
			else if(pid==2212 && pidprime==2212){
				BFy_pp[iy]-=qqprime*weight;
				BFphi_pp[iphi]-=qqprime*weight;
				DENOMy_pp[iy]+=1.0;
				DENOMphi_pp[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_pp[iqinv]-=qqprime*weight;
					DENOMqinv_pp[iqinv]+=1.0;
				}
			}
		}
	}
}

void CBF::WriteResults(int run_number){
	char filename[100];
	sprintf(filename,"results/bf%d_y.dat",run_number);
	FILE *fptr=fopen(filename,"w");
	for(int iy=0;iy<NYBINS;iy++){
		fprintf(fptr,"%7.3f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",(iy+0.5)*DELY,
		BFy_pipi[iy]/DENOMy_pipi[iy],BFy_piK[iy]/DENOMy_piK[iy],BFy_pip[iy]/DENOMy_pip[iy],
		BFy_KK[iy]/DENOMy_KK[iy],BFy_Kp[iy]/DENOMy_Kp[iy],BFy_pp[iy]/DENOMy_pp[iy]);
	}
	fclose(fptr);
	sprintf(filename,"results/bf%d_phi.dat",run_number);
	fptr=fopen(filename,"w");
	for(int iphi=0;iphi<NPHIBINS;iphi++){
		fprintf(fptr,"%7.3f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",(iphi+0.5)*DELPHI*180.0/PI,
		BFphi_pipi[iphi]/DENOMphi_pipi[iphi],BFphi_piK[iphi]/DENOMphi_piK[iphi],BFphi_pip[iphi]/DENOMphi_pip[iphi],
		BFphi_KK[iphi]/DENOMphi_KK[iphi],BFphi_Kp[iphi]/DENOMphi_Kp[iphi],BFphi_pp[iphi]/DENOMphi_pp[iphi]);
	}
	fclose(fptr);
	sprintf(filename,"results/bf%d_qinv.dat",run_number);
	fptr=fopen(filename,"w");
	for(int iqinv=0;iqinv<NQINVBINS;iqinv++){
		fprintf(fptr,"%7.3f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",(iqinv+0.5)*DELQINV,
		BFqinv_pipi[iqinv]/DENOMqinv_pipi[iqinv],BFqinv_piK[iqinv]/DENOMqinv_piK[iqinv],BFqinv_pip[iqinv]/DENOMqinv_pip[iqinv],
		BFqinv_KK[iqinv]/DENOMqinv_KK[iqinv],BFqinv_Kp[iqinv]/DENOMqinv_Kp[iqinv],BFqinv_pp[iqinv]/DENOMqinv_pp[iqinv]);
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
