#include "balhbt.h"
//using namespace std;

void CBF::Evaluate(vector<CHBTPart *> &partvec,vector<vector<CHBTPart *>> &productvec,
vector<CHBTPart *> &partprimevec,vector<vector<CHBTPart *>> &productprimevec,
vector<vector<double>> &bfnorm,unsigned int id0,unsigned int id1,unsigned int id0prime,unsigned int id1prime){
	double weight,psisquared00,psisquared01,psisquared10,psisquared11,eff,effprime;
	vector<vector<double>> psisquared(2, vector<double>(2));
	unsigned int i,iprime,iprod,iprodprime;
	int Q,Qprime,i0,i1,i0prime,i1prime;
	CHBTPart *part,*partprime;
	if(CHEAPPSISQUARED){
		psisquared[0][0]=CheapPsiSquared(partvec[0],partprimevec[0])-1.0;
		psisquared[0][1]=CheapPsiSquared(partvec[0],partprimevec[1])-1.0;
		psisquared[1][0]=CheapPsiSquared(partvec[1],partprimevec[0])-1.0;
		psisquared[1][1]=CheapPsiSquared(partvec[1],partprimevec[1])-1.0;
	}
	else{
		psisquared[0][0]=hbtcalc->GetPsiSquared(partvec[0],partprimevec[0],id0,id0prime)-1.0;
		psisquared[0][1]=hbtcalc->GetPsiSquared(partvec[0],partprimevec[1],id0,id1prime)-1.0;
		psisquared[1][0]=hbtcalc->GetPsiSquared(partvec[1],partprimevec[0],id1,id0prime)-1.0;
		psisquared[1][1]=hbtcalc->GetPsiSquared(partvec[1],partprimevec[1],id1,id1prime)-1.0;
	}
	
	for(i=0;i<2;i++){
		for(iprime=0;iprime<2;iprime++){
			if(i==0 && iprime==0){
				psisquared00=psisquared[0][0];
				psisquared01=psisquared[0][1];
				psisquared10=psisquared[1][0];
				psisquared11=psisquared[1][1];
				i0=id0;
				i1=id1;
				i0prime=id0prime;
				i1prime=id1prime;
			}
			else if(i==0 && iprime==1){
				psisquared00=psisquared[0][1];
				psisquared01=psisquared[0][0];
				psisquared10=psisquared[1][1];
				psisquared11=psisquared[1][0];
				i0=id0;
				i1=id1;
				i0prime=id1prime;
				i1prime=id0prime;
			}
			if(i==1 && iprime==0){
				psisquared00=psisquared[1][0];
				psisquared01=psisquared[1][1];
				psisquared10=psisquared[0][0];
				psisquared11=psisquared[0][1];
				i0=id1;
				i1=id0;
				i0prime=id0prime;
				i1prime=id1prime;
			}
			else{
				psisquared00=psisquared[1][1];
				psisquared01=psisquared[1][0];
				psisquared10=psisquared[0][1];
				psisquared11=psisquared[0][0];
				i0=id1;
				i1=id0;
				i0prime=id1prime;
				i1prime=id0prime;
			}
			if(UseAllWFsForCF){
				weight=psisquared00+psisquared01*bfnorm[i1prime][i0prime]+psisquared10*bfnorm[i1][i0]
					+psisquared11*bfnorm[i1][i0]*bfnorm[i1prime][i0prime];
			}
			else
				weight=psisquared00;
			Q=partvec[i]->resinfo->charge;
			Qprime=partprimevec[iprime]->resinfo->charge;
			if(Q*Qprime==1){
				NETWEIGHTsame+=weight;
			}
			if(Q*Qprime==-1){
				NETWEIGHTopp+=weight;
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
							if(abs(part->resinfo->code)==abs(partprime->resinfo->code)){
								IncrementCF(part,partprime,weight,eff*effprime);
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
	//weight=weight*efficiency;
	
	Misc::outsidelong_lcms(part->p,partprime->p,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi); // returns qout in pair frame, qout_lcms is in LCMS
	dely=fabs(dely); delphi=fabs(delphi);
	qout=fabs(qout); qside=fabs(qside); qlong=fabs(qlong); qout_lcms=fabs(qout_lcms);

	if(abs(pid)==211 && abs(pidprime)==211){
		if(pid*pidprime>0)
			cfarrays=CF_pipluspiplus;
		else{
			cfarrays=CF_pipluspiminus;
		}
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
		Ntest+=1;
		if(pid*pidprime>0)
			cfarrays=CF_pp;
		else
			cfarrays=CF_ppbar;
	}
	else{
		cfarrays=NULL;
	}
	if(cfarrays!=NULL)
		cfarrays->Increment(dely,delphi,qinv,qout_lcms,qside,qlong,weight);
}


