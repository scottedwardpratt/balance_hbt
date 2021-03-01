#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
#include "balhbt.h"

using namespace std;
int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/respars.dat");
	CResList *reslist=new CResList(&parmap);
	CResInfo::randy->reset(-1234567);
	double Tchem=150.0,taumax=100.0,strangecontent,udcontent;
	vector<vector<double>> bfnormweight;
	vector<CStableInfo *> stablevec;
	unsigned int id,id1,id2,id1prime,id2prime,NID;
	CHBTPart part1,part2,part1prime,part2prime;
	int imc,NMC=100000;
	CblastWave *bw;
	CHBTCalc hbtcalc;
	vector<CHBTPart *> partvec;
	vector<CHBTPart *> partvecprime;
	int NYBINS=20,NPHIBINS=18;
	CBF bf(NYBINS,NPHIBINS);
	
	reslist->Tf=Tchem;
	reslist->CalcEoSandChiandQdens(reslist->Tf,reslist->Pf,reslist->epsilonf,reslist->nf,reslist->densityf,
	reslist->maxweightf,reslist->chif,strangecontent,udcontent);
	reslist->chiinvf=(reslist->chif).inverse();
	reslist->FindFinalProducts(taumax);
	
	bw=new CblastWave(parmap,CResInfo::randy,reslist);
	
	GetStableInfo(reslist,taumax,stablevec,bfnormweight);
	NID=stablevec.size();
	for(id1=0;id1<NID;id1++){
		for(id2=0;id2<NID;id2++){
			printf("%10.4f ",bfnormweight[id1][id2]);
		}
		printf("\n");
	}
	
	double balweight,hbtweight;
	for(imc=0;imc<NMC;imc++){
		partvec.resize(2);
		partvecprime.resize(2);
		for(id=0;id<2;id++){
			partvec[id]=new CHBTPart();
			partvecprime[id]=new CHBTPart();
		}
		GetPart(stablevec,id1);
		GetPart(stablevec,id2);
		GetPart(stablevec,id1prime);
		GetPart(stablevec,id2prime);
		balweight=bfnormweight[id1][id2]*bfnormweight[id1prime][id2prime];
		partvec[0]->resinfo=stablevec[id1]->resinfo;
		partvec[1]->resinfo=stablevec[id2]->resinfo;
		partvecprime[0]->resinfo=stablevec[id1prime]->resinfo;
		partvecprime[1]->resinfo=stablevec[id2prime]->resinfo;
		bw->GetXP(partvec);
		bw->GetXP(partvecprime);
		
		hbtweight=hbtcalc.GetPsiSquared(partvec[0],partvecprime[0])
			*hbtcalc.GetPsiSquared(partvec[0],partvecprime[1])
				*hbtcalc.GetPsiSquared(partvec[1],partvecprime[0])
					*hbtcalc.GetPsiSquared(partvec[1],partvecprime[1]);
		
		bf.Increment(partvec,partvecprime,balweight*hbtweight);		
		
		for(id=0;id<partvec.size();id++)
			delete partvec[id];
		for(id=0;id<partvecprime.size();id++)
			delete partvecprime[id];
		partvec.clear();
		partvecprime.clear();
	}
	
	return 0;
}