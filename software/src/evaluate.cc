#include "balhbt.h"
//using namespace std;

void CBF::Evaluate(vector<CHBTPart *> &partvec,vector<vector<CHBTPart *>> &productvec,
vector<CHBTPart *> &partprimevec,vector<vector<CHBTPart *>> &productprimevec,int id0,int id1,int id0prime,int id1prime,double balweight,double balweightprime){
	double weight,psisquared00,psisquared01,psisquared10,psisquared11,eff,effprime;
	vector<vector<double>> psisquared(2, vector<double>(2));
	unsigned int i,iprime,iprod,iprodprime;
	int Q,Qprime;
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
			}
			else if(i==0 && iprime==1){
				psisquared00=psisquared[0][1];
				psisquared01=psisquared[0][0];
				psisquared10=psisquared[1][1];
				psisquared11=psisquared[1][0];
			}
			if(i==1 && iprime==0){
				psisquared00=psisquared[1][0];
				psisquared01=psisquared[1][1];
				psisquared10=psisquared[0][0];
				psisquared11=psisquared[0][1];
			}
			else{
				psisquared00=psisquared[1][1];
				psisquared01=psisquared[1][0];
				psisquared10=psisquared[0][1];
				psisquared11=psisquared[0][0];
			}
			if(UseAllWFsForCF){
				weight=psisquared00+psisquared01*balweightprime+psisquared10*balweight
					+psisquared11*balweight*balweightprime;
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
	if(cfarrays!=NULL){
		cfarrays->Increment(dely,delphi,qinv,qout_lcms,qside,qlong,weight);
	}
}

void CBF::IncrementBFs(vector<CHBTPart *> partvec,vector<vector<CHBTPart *>> productvec,int id1,int id2,double balweight){
	CHBTPart *part1,*part2;
	int pid1,pid2;
	CCF_Arrays *cfarrays;
	double eff,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi;
	unsigned int iprod1,iprod2;
	for(iprod1=0;iprod1<productvec[0].size();iprod1++){
		part1=productvec[0][iprod1];
		pid1=part1->resinfo->code;
		if(acceptancebal->acceptance(part1,eff)){
			if(abs(part1->resinfo->code)==211)
				picount+=1;
			if(abs(part1->resinfo->code)==321)
				Kcount+=1;
			if(abs(part1->resinfo->code)==2212)
				pcount+=1;
			for(iprod2=0;iprod2<productvec[1].size();iprod2++){
				part2=productvec[1][iprod2];
				pid2=part2->resinfo->code;
				if(abs(pid1)==211 && abs(pid2)==211){
					cfarrays=CF_pipluspiplus;
				}
				else if(abs(pid1)==211 && abs(pid2)==321){
					cfarrays=CF_piplusKplus;
				}
				else if(abs(pid1)==211 && abs(pid2)==2212){
					cfarrays=CF_piplusp;
				}
				else if(abs(pid1)==321 && abs(pid2)==321){
					cfarrays=CF_KplusKplus;
				}
				else if(abs(pid1)==321 && abs(pid2)==2212){
					cfarrays=CF_Kplusp;
				}
				else if(abs(pid1)==2212 && abs(pid2)==2212){
					cfarrays=CF_pp;
				}
				else{
					cfarrays=NULL;
				}
				if(cfarrays!=NULL){
					Misc::outsidelong_lcms(part1->p,part2->p,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi);
					cfarrays->Increment(dely,delphi,qinv,qout_lcms,qside,qlong,-balweight*part1->resinfo->charge*part2->resinfo->charge);

				}
			}
			for(iprod2=0;iprod2<productvec[0].size();iprod2++){
				if(iprod2!=iprod1){
					part2=productvec[0][iprod2];
					pid2=part2->resinfo->code;
					if(abs(pid1)==211 && abs(pid2)==211){
						cfarrays=CF_pipluspiplus;
					}
					else if(abs(pid1)==211 && abs(pid2)==321){
						cfarrays=CF_piplusKplus;
					}
					else if(abs(pid1)==211 && abs(pid2)==2212){
						cfarrays=CF_piplusp;
					}
					else if(abs(pid1)==321 && abs(pid2)==321){
						cfarrays=CF_KplusKplus;
					}
					else if(abs(pid1)==321 && abs(pid2)==2212){
						cfarrays=CF_Kplusp;
					}
					else if(abs(pid1)==2212 && abs(pid2)==2212){
						cfarrays=CF_pp;
					}
					else{
						cfarrays=NULL;
					}
					if(cfarrays!=NULL){
						Misc::outsidelong_lcms(part1->p,part2->p,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi);
						cfarrays->Increment(dely,delphi,qinv,qout_lcms,qside,qlong,
						-balweight*part1->resinfo->charge*part2->resinfo->charge);
					}
				}
			}
		}
	}
}

void CBF::IncrementBFCFDenoms(vector<CHBTPart *> partvec,vector<vector<CHBTPart *>> productvec,int id1,int id2,double balweight){
	balweight=1.0;
	CHBTPart *part1,*part2;
	int pid1,pid2;
	CCF_Arrays *cfarrays;
	double eff,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi;
	unsigned int iprod1,iprod2;
	for(iprod1=0;iprod1<productvec[0].size();iprod1++){
		part1=productvec[0][iprod1];
		pid1=part1->resinfo->code;
		if(acceptancebal->acceptance(part1,eff)){
			if(abs(part1->resinfo->code)==211)
				picount+=1;
			if(abs(part1->resinfo->code)==321)
				Kcount+=1;
			if(abs(part1->resinfo->code)==2212)
				pcount+=1;
			for(iprod2=0;iprod2<productvec[1].size();iprod2++){
				part2=productvec[1][iprod2];
				pid2=part2->resinfo->code;
				if(abs(pid1)==211 && abs(pid2)==211){
					cfarrays=CF_pipluspiplus;
				}
				else if(abs(pid1)==211 && abs(pid2)==321){
					cfarrays=CF_piplusKplus;
				}
				else if(abs(pid1)==211 && abs(pid2)==2212){
					cfarrays=CF_piplusp;
				}
				else if(abs(pid1)==321 && abs(pid2)==321){
					cfarrays=CF_KplusKplus;
				}
				else if(abs(pid1)==321 && abs(pid2)==2212){
					cfarrays=CF_Kplusp;
				}
				else if(abs(pid1)==2212 && abs(pid2)==2212){
					cfarrays=CF_pp;
				}
				else{
					cfarrays=NULL;
				}
				if(cfarrays!=NULL){
					Misc::outsidelong_lcms(part1->p,part2->p,qinv,qout,qout_lcms,qside,qlong,deleta,dely,delphi);
					cfarrays->Increment(dely,delphi,qinv,qout_lcms,qside,qlong,balweight);

				}
			}
		}
	}
}

