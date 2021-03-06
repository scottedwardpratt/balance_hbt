#include "balhbt.h"
using namespace std;

void CBF::Evaluate(vector<CHBTPart *> &partvec,vector<vector<CHBTPart *>> &productvec,
vector<CHBTPart *> &partprimevec,vector<vector<CHBTPart *>> &productprimevec,double balweight,double balweightprime,
int id1,int id2,int id1prime,int id2prime){
	double weight,psisquared00,psisquared01,psisquared10,psisquared11,eff,effprime,cfweight;
	unsigned int i,iprime,iprod,iprodprime;
	CHBTPart *part,*partprime;

	psisquared00=hbtcalc->GetPsiSquared(partvec[0],partprimevec[0],id1,id1prime)-1.0;
	psisquared01=hbtcalc->GetPsiSquared(partvec[0],partprimevec[1],id1,id2prime)-1.0;
	psisquared10=hbtcalc->GetPsiSquared(partvec[1],partprimevec[0],id2,id1prime)-1.0;
	psisquared11=hbtcalc->GetPsiSquared(partvec[1],partprimevec[1],id2,id2prime)-1.0;
	/*
	psisquared00=CheapPsiSquared(partvec[0],partprimevec[0])-1.0;
	psisquared01=CheapPsiSquared(partvec[0],partprimevec[1])-1.0;
	psisquared10=CheapPsiSquared(partvec[1],partprimevec[0])-1.0;
	psisquared11=CheapPsiSquared(partvec[1],partprimevec[1])-1.0;
	*/
	
	
		
	for(i=0;i<2;i++){
		for(iprime=0;iprime<2;iprime++){
			if(i==0 && iprime==0){
				weight=psisquared00+psisquared01*balweightprime+psisquared10*balweight+psisquared11*balweight*balweightprime;
				netweight+=weight*partvec[0]->resinfo->charge*partprimevec[0]->resinfo->charge;
			}
			else if(i==0 && iprime==1){
				weight=psisquared01+psisquared00*balweightprime+psisquared11*balweight+psisquared10*balweight*balweightprime;
				netweight+=weight*partvec[0]->resinfo->charge*partprimevec[1]->resinfo->charge;
			}
			else if(i==1 && iprime==0){
				weight=psisquared10+psisquared11*balweightprime+psisquared00*balweight+psisquared01*balweight*balweightprime;
				netweight+=weight*partvec[1]->resinfo->charge*partprimevec[0]->resinfo->charge;
			}
			else{
				weight=psisquared11+psisquared10*balweightprime+psisquared01*balweight+psisquared00*balweight*balweightprime;
				netweight+=weight*partvec[1]->resinfo->charge*partprimevec[1]->resinfo->charge;
			}
			for(iprod=0;iprod<productvec[i].size();iprod++){
				part=productvec[i][iprod];
				if(acceptancebal->acceptance(part,eff)){
					if(abs(part->resinfo->code)==211)
						picount+=1;
					if(abs(part->resinfo->code)==321)
						Kcount+=1;
					if(abs(part->resinfo->code)==2212)
						pcount+=1;
					for(iprodprime=0;iprodprime<productprimevec[iprime].size();iprodprime++){
						partprime=productprimevec[iprime][iprodprime];
						if(acceptancebal->acceptance(partprime,effprime)){
							Increment(part,partprime,weight,eff*effprime);
							if(abs(part->resinfo->code)==abs(partprime->resinfo->code)){
								cfweight=1.0;
								if(i==0 && iprime==0)
									cfweight=psisquared00;
								else if(i==0 && iprime==1)
									cfweight=psisquared01;
								else if(i==1 && iprime==0)
									cfweight=psisquared10;
								else if(i==1 && iprime==1)
									cfweight=psisquared11;
								IncrementCF(part,partprime,cfweight,eff*effprime);
							}
						}
					}
				}
			}
		}
	}
}

void CBF::IncrementCF(CHBTPart *part,CHBTPart *partprime,double weight,double efficiency){
	int pid,pidprime,iqinv;
	double qinv;
	pid=part->resinfo->code;
	pidprime=partprime->resinfo->code;
	if(abs(pid)!=abs(pidprime)){
		printf("In CBF::IncrementCF, |pid| != |pidprime|, %d != %d\n",abs(pid),abs(pidprime));
		exit(1);
	}
	qinv=Getqinv(part,partprime);
	iqinv=lrint(floor(qinv)/DELQINV);
	if(iqinv<NQINVBINS){
		if(abs(pid)==211){
			if(pid*pidprime<0){
				CFqinv_pipluspiminus[iqinv]+=weight;
				CF_DENOMqinv_pipluspiminus[iqinv]+=1.0;
			}
			else{
				CFqinv_pipluspiplus[iqinv]+=weight;
				CF_DENOMqinv_pipluspiplus[iqinv]+=1.0;
			}
		}
		if(abs(pid)==321){
			if(pid*pidprime<0){
				CFqinv_KplusKminus[iqinv]+=weight;
				CF_DENOMqinv_KplusKminus[iqinv]+=1.0;
			}
			else{
				CFqinv_KplusKplus[iqinv]+=weight;
				CF_DENOMqinv_KplusKplus[iqinv]+=1.0;
			}
		}
		if(abs(pid)==2212){
			if(pid*pidprime<0){
				CFqinv_ppbar[iqinv]+=weight;
				CF_DENOMqinv_ppbar[iqinv]+=1.0;
			}
			else{
				CFqinv_pp[iqinv]+=weight;
				CF_DENOMqinv_pp[iqinv]+=1.0;
			}
		}
	}
}

void CBF::Increment(CHBTPart *part,CHBTPart *partprime,double weight,double efficiency){
	int iy,iphi,iqinv,pid,pidprime;
	double y,yprime,phi,phiprime,dy,dphi,qqprime,qinv;
	pid=part->resinfo->code;
	pidprime=partprime->resinfo->code;
	qqprime=part->resinfo->charge*partprime->resinfo->charge;
	pid=abs(pid);
	pidprime=abs(pidprime);
	if((pid==211 || pid==321 || pid==2212) && (pidprime==211 || pidprime==321 || pidprime==2212)){
		//qqprime=1.0;
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
					BFqinv_pipi[iqinv]-=qqprime*weight*efficiency;
					DENOMqinv_pipi[iqinv]+=efficiency;
				}
			}
			else if((pid==211 && pidprime==321) || (pid==321 && pidprime==211)){
				BFy_piK[iy]-=qqprime*weight;
				BFphi_piK[iphi]-=qqprime*weight;
				DENOMy_piK[iy]+=1.0;
				DENOMphi_piK[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_piK[iqinv]-=qqprime*weight*efficiency;
					DENOMqinv_piK[iqinv]+=efficiency;
				}
			}
			else if((pid==211 && pidprime==2212) || (pid==2212 && pidprime==211)){
				BFy_pip[iy]-=qqprime*weight;
				BFphi_pip[iphi]-=qqprime*weight;
				DENOMy_pip[iy]+=1.0;
				DENOMphi_pip[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_pip[iqinv]-=qqprime*weight*efficiency;
					DENOMqinv_pip[iqinv]+=efficiency;
				}
			}
			else if(pid==321 && pidprime==321){
				BFy_KK[iy]-=qqprime*weight;
				BFphi_KK[iphi]-=qqprime*weight;
				DENOMy_KK[iy]+=1.0;
				DENOMphi_KK[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_KK[iqinv]-=qqprime*weight*efficiency;
					DENOMqinv_KK[iqinv]+=efficiency;
				}
			}
			else if((pid==321 && pidprime==2212) || (pid==2212 && pidprime==321)){
				BFy_Kp[iy]-=qqprime*weight;
				BFphi_Kp[iphi]-=qqprime*weight;
				DENOMy_Kp[iy]+=1.0;
				DENOMphi_Kp[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_Kp[iqinv]-=qqprime*weight*efficiency;
					DENOMqinv_Kp[iqinv]+=efficiency;
				}
			}
			else if(pid==2212 && pidprime==2212){
				BFy_pp[iy]-=qqprime*weight;
				BFphi_pp[iphi]-=qqprime*weight;
				DENOMy_pp[iy]+=1.0;
				DENOMphi_pp[iphi]+=1.0;
				if(iqinv<NQINVBINS){
					BFqinv_pp[iqinv]-=qqprime*weight*efficiency;
					DENOMqinv_pp[iqinv]+=efficiency;
				}
			}
		}
	}
}

