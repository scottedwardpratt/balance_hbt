#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
//#include "b3d.h"

class CStableInfo{
public:
	double density;
	CResInfo *resinfo;
	CStableInfo(CResInfo *resinfoset){
		resinfo=resinfoset;
		density=0.0;
	}
};

using namespace std;
int main(){
	int id,id1,id2,NID,ires,pid1,pid2;
	double strangecontent,udcontent,weight;
	CparameterMap parmap;
	CStableInfo *stableinfo;
	parmap.ReadParsFromFile("parameters/respars.dat");
	CResList *reslist=new CResList(&parmap);
	double Tchem=150.0,taumax=100.0;
	vector<vector<double>> bfnorm;
	vector<int> pid;
	vector <CStableInfo *> stablevec;
	reslist->Tf=Tchem;
	reslist->CalcEoSandChiandQdens(reslist->Tf,reslist->Pf,reslist->epsilonf,reslist->nf,reslist->densityf,
	reslist->maxweightf,reslist->chif,strangecontent,udcontent);
	reslist->chiinvf=(reslist->chif).inverse();
	reslist->FindFinalProducts(taumax);
	
	CResInfoMap::iterator rpos;
	CResInfo *resinfo,resinfoprime;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		id=0;
		if(!resinfo->decay || (HBARC/resinfo->width) > taumax){
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
	printf("NID=%d\n",NID);
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
		printf("%5d: density=%g\n",pid1,stablevec[id1]->density);
		for(id2=0;id2<NID;id2++){
			pid2=stablevec[id2]->resinfo->code;
			bfnorm[id1][id2]=reslist->CalcBalanceNorm(pid1,pid2,taumax);
			if(id1>=id2){
				printf("bfnorm[%d][%d]*dens=%g =? %g\n",id1,id2,
				bfnorm[id1][id2]*stablevec[id2]->density,bfnorm[id2][id1]*stablevec[id1]->density);
			}
		}
		printf("--------------------------------------\n");
	}
	
	double netcharge,netbaryon,netstrange;
	for(id2=0;id2<NID;id2++){
		netcharge=netbaryon=netstrange=0.0;
		for(id1=0;id1<NID;id1++){
			netcharge+=bfnorm[id1][id2]*stablevec[id1]->resinfo->charge;
			netbaryon+=bfnorm[id1][id2]*stablevec[id1]->resinfo->baryon;
			netstrange+=bfnorm[id1][id2]*stablevec[id1]->resinfo->strange;
		}
		printf("%5d: netQ=%8.5f, netB=%8.5f, netS=%8.5f\n",stablevec[id2]->resinfo->code,netcharge,netbaryon,netstrange);
	}
	return 0;
}