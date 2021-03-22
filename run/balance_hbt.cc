#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
#include "balhbt.h"

using namespace std;
int main(int argc,char *argv[]){
	int run_number;
	if (argc != 2) {
		printf("Usage: balance_hbt run_number\n");
		exit(-1);
	}
	else{
		run_number=atoi(argv[1]);
	}
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/respars.txt");
	parmap.ReadParsFromFile("parameters/bfpars.txt");
	CBalHBT *balhbt=new CBalHBT(&parmap,run_number);

	double Tchem=150.0,taumax=100.0,strangecontent,udcontent,balweight;
	vector<vector<double>> bfnorm;
	vector<CStableInfo *> stablevec;
	unsigned int id,id1,id2,id1prime,id2prime,NID,i,iprod;
	CHBTPart part1,part2,part1prime,part2prime;
	double nhadron0;
	long long int imc,NMC=10000;
	printf("Enter NMC: ");
	scanf("%lld",&NMC);

	balhbt->reslist->Tf=Tchem;
	balhbt->reslist->CalcEoSandChiandQdens(balhbt->reslist->Tf,balhbt->reslist->Pf,balhbt->reslist->epsilonf,balhbt->reslist->nf,balhbt->reslist->densityf,
	balhbt->reslist->maxweightf,balhbt->reslist->chif,strangecontent,udcontent);
	balhbt->reslist->FindFinalProducts(taumax);
	balhbt->bw=new CblastWave(&parmap,balhbt->randy,balhbt->reslist);
	balhbt->GetStableInfo(balhbt->reslist,taumax,stablevec,bfnorm);
	NID=stablevec.size();
	balhbt->InitHBT(stablevec,"parameters/hbtpars.txt");
	nhadron0=balhbt->reslist->nf;
	
	for(id1=0;id1<NID;id1++){
		for(id2=0;id2<NID;id2++){
			printf("%10.4f ",bfnorm[id1][id2]);
		}
		printf("\n");
	}
	vector<CHBTPart *> partvec(2);
	vector<CHBTPart *> partprimevec(2);
	for(id=0;id<2;id++){
		partvec[id]=new CHBTPart();
		partprimevec[id]=new CHBTPart();
	}
	vector<vector<CHBTPart *>> productvec(2);
	vector<vector<CHBTPart *>> productprimevec(2);
	
	for(imc=0;imc<NMC;imc++){
		if(imc%(NMC/10)==0)
			printf("finished %ld percent\n",lrint(100.0*imc/double(NMC)));
		balhbt->GetPart(stablevec,id1);
		balhbt->GetPart(stablevec,id2);
		balhbt->GetPart(stablevec,id1prime);
		balhbt->GetPart(stablevec,id2prime);
		balweight=bfnorm[id1][id2]*bfnorm[id1prime][id2prime]*nhadron0*nhadron0;
		partvec[0]->resinfo=stablevec[id1]->resinfo;
		partvec[1]->resinfo=stablevec[id2]->resinfo;
		partprimevec[0]->resinfo=stablevec[id1prime]->resinfo;
		partprimevec[1]->resinfo=stablevec[id2prime]->resinfo;
		
		balhbt->bw->GetXP(partvec);
		balhbt->bw->GetXP(partprimevec);
		
		balhbt->GetDecayProducts(partvec[0],productvec[0]);
		balhbt->GetDecayProducts(partvec[1],productvec[1]);
		balhbt->GetDecayProducts(partprimevec[0],productprimevec[0]);
		balhbt->GetDecayProducts(partprimevec[1],productprimevec[1]);

		balhbt->bf->Evaluate(partvec,productvec,partprimevec,productprimevec,balweight);
		
		for(i=0;i<2;i++){
			for(iprod=0;iprod<productvec[i].size();iprod++)
				delete productvec[i][iprod];
			productvec[i].clear();
			for(iprod=0;iprod<productprimevec[i].size();iprod++)
				delete productprimevec[i][iprod];
			productprimevec[i].clear();	
		}
		balhbt->bf->WriteResults(run_number);
	}

	return 0;
}