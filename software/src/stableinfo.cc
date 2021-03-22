#include "balhbt.h"
using namespace std;
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
			if(resinfo->code!=22 && resinfo->code>0 && (resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0)){
				stableinfo=new CStableInfo(resinfo);
				stablevec.push_back(stableinfo);
				//if(resinfo->decay)
				printf("PID=%d, tau=%g, %s\n",resinfo->code,HBARC/resinfo->width,resinfo->name.c_str());
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
			printf("bfnorm[%d][%d]=%g, bfnorm[%d][%d]*dens=%g =? %g\n",id1,id2,bfnorm[id1][id2],id1,id2,
			bfnorm[id1][id2]*stablevec[id2]->density,bfnorm[id2][id1]*stablevec[id1]->density);
		}
		printf("--------------------------------------\n");
	}
	
	double netcharge,netbaryon,netstrange;
	CStableInfo::denstot=0.0;
	for(id2=0;id2<NID;id2++){
		netcharge=netbaryon=netstrange=0.0;
		for(id1=0;id1<NID;id1++){
			netcharge+=bfnorm[id1][id2]*stablevec[id1]->resinfo->charge;
			netbaryon+=bfnorm[id1][id2]*stablevec[id1]->resinfo->baryon;
			netstrange+=bfnorm[id1][id2]*stablevec[id1]->resinfo->strange;
			bfnorm[id1][id2]=bfnorm[id1][id2]/stablevec[id1]->density;
		}
		CStableInfo::denstot+=stablevec[id2]->density;
		printf("%5d: netQ=%8.5f, netB=%8.5f, netS=%8.5f\n",stablevec[id2]->resinfo->code,netcharge,netbaryon,netstrange);
	}

}

void CBalHBT::GetPart(vector<CStableInfo *> &stablevec,unsigned int &id){
	id=0;
	double denstarget=CStableInfo::denstot*CResInfo::randy->ran();
	double netdens=stablevec[id]->density;
	while(denstarget>netdens){
		id+=1;
		netdens+=stablevec[id]->density;
	}
}