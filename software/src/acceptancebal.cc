#include "balhbt.h"
//using namespace std;
CAcceptanceBal *CBF::acceptancebal=NULL;

CAcceptanceBal::CAcceptanceBal(CparameterMap *parmap){
	ACCEPTANCE=parmap->getS("BF_ACCEPTANCE","STAR");
	if(ACCEPTANCE=="STAR"){
		etamax=0.9;
		ptmin_pi=200.0;
		ptmin_K=200.0;
		ptmin_p=400.0;
	}
	else if(ACCEPTANCE=="ALICE"){
		etamax=0.8;
		ptmin_pi=200.0;
		ptmin_K=200.0;
		ptmin_p=200.0;
	}
	else{
		fprintf(balhbt->logfile,"ACCEPTANCE not recognized in CAcceptanceBal::CAcceptanceBal()\n");
		exit(0);
	}
}

bool CAcceptanceBal::acceptance(CHBTPart *part1,CHBTPart *part2,double &efficiency){
	double pt1,pt2,eta1,eta2;
	efficiency=0.0;
	if(abs(part1->resinfo->code)!=211 && abs(part1->resinfo->code)!=321 && abs(part1->resinfo->code)!=2212)
		return false;
	if(abs(part2->resinfo->code)!=211 && abs(part2->resinfo->code)!=321 && abs(part2->resinfo->code)!=2212)
		return false;
	eta1=atanh(part1->p[3]/part1->p[0]);
	eta2=atanh(part2->p[3]/part2->p[0]);
	if(fabs(eta1-eta2)>2.0*etamax){
		return false;
	}
	
	pt1=sqrt(part1->p[1]*part1->p[1]+part1->p[2]*part1->p[2]);
	if(abs(part1->resinfo->code)==211 && pt1<ptmin_pi)
		return false;
	if(abs(part1->resinfo->code)==321 && pt1<ptmin_K)
		return false;
	if(abs(part1->resinfo->code)==2212 && pt1<ptmin_p)
		return false;
		
	pt2=sqrt(part2->p[1]*part2->p[1]+part2->p[2]*part2->p[2]);
	if(abs(part2->resinfo->code)==211 && pt2<ptmin_pi)
		return false;
	if(abs(part2->resinfo->code)==321 && pt2<ptmin_K)
		return false;
	if(abs(part2->resinfo->code)==2212 && pt2<ptmin_p)
		return false;
	
	efficiency=1.0;
	if(ACCEPTANCE=="STAR")
		efficiency=1.0-0.5*fabs(eta1-eta2)/etamax;

	return true;	
}

bool CAcceptanceBal::acceptance(CHBTPart *part,double &efficiency){
	bool accept=false;
	int abspid=abs(part->resinfo->code);
	efficiency=0.0;
	double dca=part->GetDCA();
	if(dca>3.0)
		return false;
	if(abspid!=211 && abspid!=321 && abspid!=2212)
		return false;
	else if(fabs(part->eta)>etamax)
		return false;
	else if(abspid==211 && part->pt<ptmin_pi)
		return false;
	else if(abspid==321 && part->pt<ptmin_K)
		return false;
	else if(abspid==2212 && part->pt<ptmin_p)
		return false;
	else{
		efficiency=1.0;
		accept=true;
	}
	return accept;
}

bool CAcceptanceBal::acceptance_spectra(CHBTPart *part,double &efficiency){
	bool accept=false;
	double dca=part->GetDCA();
	if(fabs(part->y)<0.5 && fabs(dca)<0.001)
		accept=true;
	return accept;
}