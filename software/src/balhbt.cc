#include "balhbt.h"
//using namespace std;
CBF *CBalHBT::bf=NULL;

CBalHBT::CBalHBT(int run_number_set){
	run_number=run_number_set;
	parmap.ReadParsFromFile("parameters/respars.txt");
	parmap.ReadParsFromFile("parameters/bwpars.txt");
	parmap.ReadParsFromFile("parameters/bfpars.txt");
}

CBalHBT::CBalHBT(int run_number_set,double BW_T,double BW_UPERP){
	run_number=run_number_set;
	parmap.ReadParsFromFile("parameters/bwpars.txt");
	parmap.ReadParsFromFile("parameters/respars.txt");
	parmap.ReadParsFromFile("parameters/bfpars.txt");
	parmap.set("BW_T",BW_T);
	parmap.set("BW_UPERP",BW_UPERP);
}

void CBalHBT::Init(){
	double strangecontent,udcontent;
	string logfilename="logfiles/balhbt"+to_string(run_number)+".txt";
	logfile=fopen(logfilename.c_str(),"a");
	
	NMC=parmap.getI("BF_NMC",1000000);
	BARYONSONLY=parmap.getB("BF_BARYONSONLY",false); // only baryons (for pp correlations)
	STRANGEONLY=parmap.getB("BF_STRANGEONLY",false); // only strange particles (for KK correlations)
	Tchem=parmap.getD("BF_TCHEM",150.0);
	CBalHBT::bf=new CBF(this);
	reslist=new CResList(&parmap);
	randy=new CRandy(run_number-time(NULL));
	decay_nbody=new CDecay_NBody(randy);
	hbtcalc=new CHBTCalc(&parmap);
	bf->hbtcalc=hbtcalc;
	bf->randy=randy;
	bf->balhbt=this;
	CResInfo::randy=randy;
	//
	taumax=parmap.getD("BF_TAUMAX",100.0);
	reslist->Tf=Tchem;
	reslist->CalcEoSandChiandQdens(reslist->Tf,reslist->Pf,reslist->epsilonf,
	reslist->nf,reslist->densityf,
	reslist->maxweightf,reslist->chif,strangecontent,udcontent);
	reslist->FindFinalProducts(taumax);
	bw=new CblastWave(&parmap,randy,reslist);
	GetStableInfo(reslist,taumax,stablevec,bfnorm);
	fprintf(logfile,"Initialized blast wave, resonance and decay info for balhbt\n");
	fflush(logfile);
	bool SPECTRA_ONLY=parmap.getB("BF_SPECTRA_ONLY",false);
	if(!SPECTRA_ONLY)
		InitHBT(stablevec,"parameters/hbtpars.txt");
};

CBalHBT::~CBalHBT(){
	fclose(logfile);
	delete reslist;
	delete randy;
	delete bf;
	delete hbtcalc;
}

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
					symweight=0.5;
				}
				hbtcalc->wf_same[id1][id2]=new CWaveFunction_generic(parsfilename,q1q2,info1->resinfo->mass,info2->resinfo->mass,symweight);
			}
		}
	}
	hbtcalc->wf_pp=new CWaveFunction_pp_schrod(parsfilename);
	hbtcalc->wf_classical=new CWaveFunction_classical();
	fprintf(logfile,"Initialized Wave Functions for HBT\n");
	fflush(logfile);
}

void CBalHBT::CalcCFs(){
	long long unsigned int imc;
	unsigned int id,id1,id2,i,iprod,id1prime,id2prime;
	double balweight,balweightprime;
	vector<CHBTPart *> partvec(2);
	vector<CHBTPart *> partprimevec(2);
	vector<vector<CHBTPart *>> productvec(2);
	vector<vector<CHBTPart *>> productprimevec(2);
	for(id=0;id<2;id++){
		partvec[id]=new CHBTPart();
		partprimevec[id]=new CHBTPart();
	}
	for(imc=0;imc<NMC;imc++){
		GetPart(stablevec,id1);
		GetPart(stablevec,id2);
		GetPart(stablevec,id1prime);
		GetPart(stablevec,id2prime);
		if(BARYONSONLY){
			while(stablevec[id1]->resinfo->baryon==0 && stablevec[id2]->resinfo->baryon==0){
				GetPart(stablevec,id1);
				GetPart(stablevec,id2);
			} 
			while(stablevec[id1prime]->resinfo->baryon==0 && stablevec[id2prime]->resinfo->baryon==0){
				GetPart(stablevec,id1prime);
				GetPart(stablevec,id2prime);
			}
		}
		else if(STRANGEONLY){
			while(stablevec[id1]->resinfo->strange==0 && stablevec[id2]->resinfo->strange==0){
				GetPart(stablevec,id1);
				GetPart(stablevec,id2);
			}
			while(stablevec[id1prime]->resinfo->strange==0 && stablevec[id2prime]->resinfo->strange==0){
				GetPart(stablevec,id1prime);
				GetPart(stablevec,id2prime);
			}
		}
		balweight=-bfnorm[id1][id2]*0.5;
		balweightprime=-bfnorm[id1prime][id2prime]*0.5;
		partvec[0]->resinfo=stablevec[id1]->resinfo;
		partvec[1]->resinfo=stablevec[id2]->resinfo;
		partprimevec[0]->resinfo=stablevec[id1prime]->resinfo;
		partprimevec[1]->resinfo=stablevec[id2prime]->resinfo;
		if(randy->ran()<0.5){
			partvec[0]->PartAntipart();
			balweight=-balweight;
		}
		if(randy->ran()<0.5){
			partvec[1]->PartAntipart();
			balweight=-balweight;
		}
		if(randy->ran()<0.5){
			partprimevec[0]->PartAntipart();
			balweightprime=-balweightprime;
		}
		if(randy->ran()<0.5){
			partprimevec[1]->PartAntipart();
			balweightprime=-balweightprime;
		}
		
		bw->GetXP(partvec);
		bw->GetXP(partprimevec);
		
		GetDecayProducts(partvec[0],productvec[0]);
		GetDecayProducts(partvec[1],productvec[1]);
		GetDecayProducts(partprimevec[0],productprimevec[0]);
		GetDecayProducts(partprimevec[1],productprimevec[1]);
		
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
			bf->Evaluate(partvec,productvec,partprimevec,productprimevec,bfnorm,id1,id2,id1prime,id2prime);
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
			fprintf(logfile,"finished %ld percent\n",lrint(100.0*imc/double(NMC)));
			fflush(logfile);
		}
	}
	
	for(id=0;id<2;id++){
		delete partvec[id];
		delete partprimevec[id];
	}
	partvec.clear();
	partprimevec.clear();
	bf->WriteResults(run_number);
	fclose(logfile);
}