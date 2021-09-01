#include "balhbt.h"
using namespace std;

void CBF::Evaluate(vector<CHBTPart *> &partvec,vector<vector<CHBTPart *>> &productvec,
vector<CHBTPart *> &partprimevec,vector<vector<CHBTPart *>> &productprimevec,double balweight,double balweightprime,
int id1,int id2,int id1prime,int id2prime){
	double weight,cfweight=0.0,psisquared00,psisquared01,psisquared10,psisquared11,eff,effprime;
	unsigned int i,iprime,iprod,iprodprime;
	CHBTPart *part,*partprime;
	if(CHEAPPSISQUARED){
		psisquared00=CheapPsiSquared(partvec[0],partprimevec[0])-1.0;
		psisquared01=CheapPsiSquared(partvec[0],partprimevec[1])-1.0;
		psisquared10=CheapPsiSquared(partvec[1],partprimevec[0])-1.0;
		psisquared11=CheapPsiSquared(partvec[1],partprimevec[1])-1.0;
	
	}
	else{
		psisquared00=hbtcalc->GetPsiSquared(partvec[0],partprimevec[0],id1,id1prime)-1.0;
		psisquared01=hbtcalc->GetPsiSquared(partvec[0],partprimevec[1],id1,id2prime)-1.0;
		psisquared10=hbtcalc->GetPsiSquared(partvec[1],partprimevec[0],id2,id1prime)-1.0;
		psisquared11=hbtcalc->GetPsiSquared(partvec[1],partprimevec[1],id2,id2prime)-1.0;
	}
	
	for(i=0;i<2;i++){
		for(iprime=0;iprime<2;iprime++){
			if(i==0 && iprime==0){
				weight=psisquared00+psisquared01*balweightprime+psisquared10*balweight+psisquared11*balweight*balweightprime;
			}
			else if(i==0 && iprime==1){
				weight=psisquared01+psisquared00*balweightprime+psisquared11*balweight+psisquared10*balweight*balweightprime;
			}
			else if(i==1 && iprime==0){
				weight=psisquared10+psisquared11*balweightprime+psisquared00*balweight+psisquared01*balweight*balweightprime;
			}
			else{
				weight=psisquared11+psisquared10*balweightprime+psisquared01*balweight+psisquared00*balweight*balweightprime;
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
							if(UseAllWFsForCF){
								cfweight=weight;
							}
							else{
								if(i==0 && iprime==0)
									cfweight=psisquared00;
								if(i==0 && iprime==1)
									cfweight=psisquared01;
								if(i==1 && iprime==0)
									cfweight=psisquared10;
								if(i==1 && iprime==1)
									cfweight=psisquared11;
							}
							if(abs(part->resinfo->code)==abs(partprime->resinfo->code)){
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
	int pid,pidprime;
	double qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi;
	CCF_Arrays *cfarrays=NULL;
	pid=part->resinfo->code;
	pidprime=partprime->resinfo->code;
	if(abs(pid)!=abs(pidprime)){
		fprintf(balhbt->logfile,"In CBF::IncrementCF, |pid| != |pidprime|, %d != %d\n",abs(pid),abs(pidprime));
		exit(0);
	}
	weight=weight*efficiency;
	
	Misc::outsidelong_lcms(part->p,partprime->p,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi); // returns qout in pair frame, qout_lcms is in LCMS
	dely=fabs(dely); delphi=fabs(delphi);
	qout=fabs(qout); qside=fabs(qside); qlong=fabs(qlong); qout_lcms=fabs(qout_lcms);
	//printf("qout=%g, qout_lcms=%g\n",qout,qout_lcms);
	if(abs(pid)==211 && abs(pidprime)==211){
		if(pid*pidprime>0)
			cfarrays=CF_pipluspiplus;
		else
			cfarrays=CF_pipluspiminus;
	}
	else if(abs(pid)==211 && abs(pidprime)==321){
		if(pid*pidprime>0)
			cfarrays=CF_piplusKplus;
		else
			cfarrays=CF_piplusKminus;
	}
	else if(abs(pid)==211 && abs(pidprime)==2212){
		if(pid*pidprime>0)
			cfarrays=CF_piplusp;
		else
			cfarrays=CF_pipluspbar;
	}
	else if(abs(pid)==321 && abs(pidprime)==321){
		if(pid*pidprime>0)
			cfarrays=CF_KplusKplus;
		else
			cfarrays=CF_KplusKminus;
	}
	else if(abs(pid)==321 && abs(pidprime)==2212){
		if(pid*pidprime>0)
			cfarrays=CF_Kplusp;
		else
			cfarrays=CF_Kpluspbar;
	}
	else if(abs(pid)==2212 && abs(pidprime)==2212){
		if(pid*pidprime>0)
			cfarrays=CF_pp;
		else
			cfarrays=CF_ppbar;
	}
	else{
		cfarrays=NULL;
	}
	cfarrays->Increment(dely,delphi,qinv,qout_lcms,qside,qlong,weight);
}


