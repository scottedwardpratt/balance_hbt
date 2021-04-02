#include "balhbt.h"
using namespace std;
#include "balhbt.h"
using namespace std;

CBalHBT::CBalHBT(CparameterMap *parmapin,int run_number){
	parmap=parmapin;
	reslist=new CResList(parmap);
	randy=new CRandy(run_number);
	bf=new CBF(parmap);
	hbtcalc=new CHBTCalc(parmap);
	bf->hbtcalc=hbtcalc;
	bf->randy=randy;
	CResInfo::randy=randy;
};

void CBalHBT::InitHBT(vector<CStableInfo *> &stablevec,string parsfilename){
	int NID=stablevec.size();
	int id1,id2,q1,q2,q1q2;
	CStableInfo *info1,*info2;
	double symweight;
	hbtcalc->wf_same.resize(NID);
	hbtcalc->wf_opp.resize(NID);
	for(id1=0;id1<NID;id1++){
		hbtcalc->wf_same[id1].resize(NID);
		hbtcalc->wf_opp[id1].resize(NID);
		for(id2=0;id2<NID;id2++){
			hbtcalc->wf_same[id1][id2]=NULL;
			hbtcalc->wf_opp[id1][id2]=NULL;
		}
	}

	for(id1=0;id1<NID;id1++){
		info1=stablevec[id1];
		q1=info1->resinfo->charge;
		for(id2=0;id2<=id1;id2++){
			info2=stablevec[id2];
			q2=info2->resinfo->charge;
			q1q2=abs(q1*q2);
			if(id1!=id2){
				if(q1q2!=0){
					hbtcalc->wf_same[id1][id2]=new CWaveFunction_generic(parsfilename,q1q2,info1->resinfo->mass,info2->resinfo->mass,0.5);
					hbtcalc->wf_same[id2][id1]=hbtcalc->wf_same[id1][id2];
					hbtcalc->wf_opp[id1][id2]=new CWaveFunction_generic(parsfilename,-q1q2,info1->resinfo->mass,info2->resinfo->mass,0.5);
					hbtcalc->wf_opp[id2][id1]=hbtcalc->wf_opp[id1][id2];
				}
			}
			else{
				hbtcalc->wf_opp[id1][id2]=new CWaveFunction_generic(parsfilename,-q1q2,info1->resinfo->mass,info2->resinfo->mass,0.5);
				if(fabs(info1->resinfo->spin)<1.0E-3)
					symweight=1.0;
				else if(fabs(info1->resinfo->spin-0.5)<1.0E-3)
					symweight=0.25;
				else if(fabs(info1->resinfo->spin-1.5)<1.0E-3)
					symweight=0.375;
				else{
					printf("spin not recognized\n");
					symweight=0.5;
				}
				hbtcalc->wf_same[id1][id2]=new CWaveFunction_generic(parsfilename,q1q2,info1->resinfo->mass,info2->resinfo->mass,symweight);
			}
		}
	}
	hbtcalc->wf_pp=new CWaveFunction_pp_schrod(parsfilename);
}