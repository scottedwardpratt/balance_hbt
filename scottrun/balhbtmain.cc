#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
#include "balhbt.h"

using namespace std;
int main(int argc,char *argv[]){
	int run_number,success=0;
	double R,tau;
	if (argc != 4) {
		printf("Usage: balance_hbt tau R run_number\n");
		exit(-1);
	}
	else{
		tau=atof(argv[1]);
		R=atof(argv[2]);
		run_number=atoi(argv[3]);
	}
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters/respars.txt");
	parmap.ReadParsFromFile("parameters/bwpars.txt");
	parmap.ReadParsFromFile("parameters/bfpars.txt");
	CBalHBT *balhbt=new CBalHBT(&parmap,run_number);
	
	double Tchem=150.0,taumax=100.0,strangecontent,udcontent,balweight,balweightprime,nhadron0;
	vector<vector<double>> bfnorm;
	vector<CStableInfo *> stablevec;
	unsigned int id,id1,id2,id1prime,id2prime,i,iprod;
	CHBTPart part1,part2,part1prime,part2prime;
	long long int imc,NMC;
	
	balhbt->reslist->Tf=Tchem;
	balhbt->reslist->CalcEoSandChiandQdens(balhbt->reslist->Tf,balhbt->reslist->Pf,balhbt->reslist->epsilonf,
	balhbt->reslist->nf,balhbt->reslist->densityf,
	balhbt->reslist->maxweightf,balhbt->reslist->chif,strangecontent,udcontent);
	balhbt->reslist->FindFinalProducts(taumax);
	balhbt->bw=new CblastWave(&parmap,balhbt->randy,balhbt->reslist);
	balhbt->bw->tau=tau;
	balhbt->bw->Rperp=R;
	balhbt->GetStableInfo(balhbt->reslist,taumax,stablevec,bfnorm);
	balhbt->InitHBT(stablevec,"parameters/hbtpars.txt");
	NMC=balhbt->parmap->getD("BF_NMC",100000);
	nhadron0=CStableInfo::denstot;
	
	vector<CHBTPart *> partvec(2);
	vector<CHBTPart *> partprimevec(2);
	for(id=0;id<2;id++){
		partvec[id]=new CHBTPart();
		partprimevec[id]=new CHBTPart();
	}
	vector<vector<CHBTPart *>> productvec(2);
	vector<vector<CHBTPart *>> productprimevec(2);
	fprintf(balhbt->logfile,"finished initialization\n");
	fflush(balhbt->logfile);
	
	for(imc=0;imc<NMC;imc++){
		balhbt->GetPart(stablevec,id1);
		balhbt->GetPart(stablevec,id2);
		balhbt->GetPart(stablevec,id1prime);
		balhbt->GetPart(stablevec,id2prime);
		balweight=-bfnorm[id1][id2]*nhadron0*0.5;
		balweightprime=-bfnorm[id1prime][id2prime]*nhadron0*0.5;
		partvec[0]->resinfo=stablevec[id1]->resinfo;
		partvec[1]->resinfo=stablevec[id2]->resinfo;
		partprimevec[0]->resinfo=stablevec[id1prime]->resinfo;
		partprimevec[1]->resinfo=stablevec[id2prime]->resinfo;
		if(balhbt->randy->ran()<0.5){
			partvec[0]->PartAntipart();
			balweight=-balweight;
		}
		if(balhbt->randy->ran()<0.5){
			partvec[1]->PartAntipart();
			balweight=-balweight;
		}
		if(balhbt->randy->ran()<0.5){
			partprimevec[0]->PartAntipart();
			balweightprime=-balweightprime;
		}
		if(balhbt->randy->ran()<0.5){
			partprimevec[1]->PartAntipart();
			balweightprime=-balweightprime;
		}
		balhbt->bw->GetXP(partvec);
		balhbt->bw->GetXP(partprimevec);
		
		balhbt->GetDecayProducts(partvec[0],productvec[0]);
		balhbt->GetDecayProducts(partvec[1],productvec[1]);
		balhbt->GetDecayProducts(partprimevec[0],productprimevec[0]);
		balhbt->GetDecayProducts(partprimevec[1],productprimevec[1]);
		
		for(int anti=0;anti<2;anti++){
			if(anti==1){
				partvec[0]->PartAntipart();
				for(unsigned long int iprod=0;iprod<productvec[0].size();iprod++){
					productvec[0][iprod]->PartAntipart();
				}
				partvec[1]->PartAntipart();
				for(unsigned long int iprod=0;iprod<productvec[1].size();iprod++){
					productvec[1][iprod]->PartAntipart();
				}
			}
			balhbt->bf->Evaluate(partvec,productvec,partprimevec,productprimevec,balweight,balweightprime,id1,id2,id1prime,id2prime);
		}
		
		for(i=0;i<2;i++){
			for(iprod=0;iprod<productvec[i].size();iprod++)
				delete productvec[i][iprod];
			productvec[i].clear();
			for(iprod=0;iprod<productprimevec[i].size();iprod++)
				delete productprimevec[i][iprod];
			productprimevec[i].clear();
		}
		if((imc+1)%(NMC/10)==0){
			fprintf(balhbt->logfile,"finished %ld percent\n",lrint(100.0*imc/double(NMC)));
			fflush(balhbt->logfile);
		}
	}
	balhbt->bf->WriteResults(run_number);
	fclose(balhbt->logfile);
	success=1;
	printf("%d\n",success);
	return 0;
}