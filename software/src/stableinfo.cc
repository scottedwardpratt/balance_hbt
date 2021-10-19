#include "balhbt.h"
//using namespace std;
double CStableInfo::denstot=0.0;

void CBalHBT::GetStableInfo(CResList *reslist,double taumax,vector<CStableInfo *> &stablevec,vector<vector<double>> &bfnorm){
	int id,id1,id2,NID,ires,pid1,pid2;
	double weight;
	CStableInfo *stableinfo;
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		id=0;
		if(!resinfo->decay || (HBARC/resinfo->width) > 100.0){
			//if(resinfo->code!=22 && resinfo->code>0 && (resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0)){
			if(resinfo->code!=22 && resinfo->code>0){
				stableinfo=new CStableInfo(resinfo);
				stablevec.push_back(stableinfo);
			}
			id+=1;
		}
	}
	NID=stablevec.size();
	bfnorm.resize(NID);
	for(id1=0;id1<NID;id1++){
		bfnorm[id1].resize(NID);
		for(id2=0;id2<NID;id2++)
			bfnorm[id1][id2]=0.0;
	}
		
	for(id1=0;id1<NID;id1++){
		pid1=stablevec[id1]->resinfo->code;
		//stablevec[id1]->density=reslist->densityf[stablevec[id1]->resinfo->ires];
		stablevec[id1]->density=0.0;
		for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
			resinfo=rpos->second;
			ires=resinfo->ires;
			if(resinfo->FindContent(pid1,1.0,taumax,weight)){
				stablevec[id1]->density+=weight*reslist->densityf[ires];
			}
		}
	}
	for(id1=0;id1<NID;id1++){
		pid1=stablevec[id1]->resinfo->code;
		for(id2=0;id2<NID;id2++){
			pid2=stablevec[id2]->resinfo->code;
			bfnorm[id1][id2]=reslist->CalcBalanceNorm(pid1,pid2,taumax);
		}
	}
	
	CStableInfo::denstot=0.0;
	for(id2=0;id2<NID;id2++){
		resinfo=stablevec[id2]->resinfo;
		if(resinfo->code>0 && resinfo->code!=22){
			if(resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0)
				CStableInfo::denstot+=2.0*stablevec[id2]->density;
			else
				CStableInfo::denstot+=stablevec[id2]->density;
		}
		/*
		double netcharge,netbaryon,netstrange;
		netcharge=netbaryon=netstrange=0.0;
		for(id1=0;id1<NID;id1++){
		netcharge+=bfnorm[id1][id2]*stablevec[id1]->resinfo->charge;
		netbaryon+=bfnorm[id1][id2]*stablevec[id1]->resinfo->baryon;
		netstrange+=bfnorm[id1][id2]*stablevec[id1]->resinfo->strange;
		}
		printf("netcharge=%8.6f=?%d, netbaryon=%8.6f=?%d, netstrange=%8.6f=?%d\n",netcharge,stablevec[id2]->resinfo->charge,
		netbaryon,stablevec[id2]->resinfo->baryon,netstrange,stablevec[id2]->resinfo->strange);*/
	}
	for(id2=0;id2<NID;id2++){
		for(id1=0;id1<NID;id1++){
			resinfo=stablevec[id2]->resinfo;
			if(resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0)
				bfnorm[id1][id2]=bfnorm[id1][id2]*CStableInfo::denstot/(2.0*stablevec[id1]->density);
			else
				bfnorm[id1][id2]=0.0;
		}
	}
}

void CBalHBT::GetPart(vector<CStableInfo *> &stablevec,unsigned int &id){
	id=0;
	CResInfo *resinfo;
	double denstarget=CStableInfo::denstot*CResInfo::randy->ran();
	if(denstarget>CStableInfo::denstot){
		fprintf(logfile,"CStableInfo::denstot=%g, denstarget=%g\n",CStableInfo::denstot,denstarget);
		exit(1);
	}
	double netdens=0.0;
	resinfo=stablevec[id]->resinfo;
	if(resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0)
		netdens+=2.0*stablevec[id]->density;
	else
		netdens+=stablevec[id]->density;
	while(denstarget>netdens){
		id+=1;
		resinfo=stablevec[id]->resinfo;
		if(resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0)
			netdens+=2.0*stablevec[id]->density;
		else
			netdens+=stablevec[id]->density;
		if(netdens>CStableInfo::denstot){
			fprintf(logfile,"YIKES: id=%d, NID=%ld\n",id,stablevec.size());
			fprintf(logfile,"id=%d: pid=%d, density=%g, CStableInfo::denstot=%g, netdens=%g\n",id,resinfo->code,stablevec[id]->density,CStableInfo::denstot,netdens);
			exit(1);
		}
	}
	
}