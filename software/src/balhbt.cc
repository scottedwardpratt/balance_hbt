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
	int id1,id2,q1,q2;
	CStableInfo *info1,*info2;
	double symweight;
	hbtcalc->wf.resize(NID);
	for(id1=0;id1<NID;id1++)
		hbtcalc->wf[id1].resize(NID);
	for(id1=0;id1<NID;id1++){
		info1=stablevec[id1];
		q1=info1->resinfo->charge;
		for(id2=0;id2<id1;id2++){
			info2=stablevec[id2];
			q2=info2->resinfo->charge;
			if(id1!=id2){
				symweight=0.5;
			}
			else{
				if(abs(info1->resinfo->code)==2212){
					hbtcalc->wf[id1][id2]=new CWaveFunction_pp_schrod(parsfilename);
				}
				else if(info1->resinfo->spin==0)
					symweight=1.0;
				else if(fabs(info1->resinfo->spin-0.5)<1.0E-5)
					symweight=0.25;
				else
					symweight=0.5;
			}
			hbtcalc->wf[id1][id2]=new CWaveFunction_generic(parsfilename,q1*q2,info1->resinfo->mass,info2->resinfo->mass,symweight);	
		}
		if(id1!=id2)
			hbtcalc->wf[id2][id1]=hbtcalc->wf[id1][id2];
	}
}