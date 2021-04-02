#include "balhbt.h"
using namespace std;

CHBTCalc::CHBTCalc(CparameterMap *parmapin){
	parmap=parmapin;
}

double CHBTCalc::GetPsiSquared(CHBTPart *part,CHBTPart *partprime,int id,int idprime){
	double q2,r2,qdotr,P2,Pdotq,Pdotr,psisquared;
	double qmag,rmag,ctheta;
	int alpha,pid,pidprime;
	FourVector q,r,P;
	pid=part->resinfo->code;
	pidprime=partprime->resinfo->code;
	int q1q2=part->resinfo->charge*partprime->resinfo->charge;
	q2=qdotr=r2=P2=0.0;
	double g[4]={1.0,-1.0,-1.0,-1.0};
	q2=qdotr=P2=Pdotr=Pdotq=0.0;
	for(alpha=0;alpha<4;alpha++){
		q[alpha]=0.5*(part->p[alpha]-partprime->p[alpha]);
		r[alpha]=part->x[alpha]-partprime->x[alpha];
		P[alpha]=part->p[alpha]+partprime->p[alpha];
		r2+=g[alpha]*r[alpha]*r[alpha];
		q2+=g[alpha]*q[alpha]*q[alpha];
		qdotr+=g[alpha]*q[alpha]*r[alpha];
		P2+=g[alpha]*P[alpha]*P[alpha];
		Pdotr+=g[alpha]*P[alpha]*r[alpha];
		Pdotq+=g[alpha]*P[alpha]*q[alpha];
	}
	q2=q2-Pdotq*Pdotq/P2;
	r2=r2-Pdotr*Pdotr/P2;
	qdotr=qdotr-(Pdotr*Pdotq)/P2;
	qmag=sqrt(-q2);
	rmag=sqrt(-r2);
	ctheta=-qdotr/(qmag*rmag);
	if(q1q2>0){
		if(abs(pid)==2212 && pid==pidprime && qmag<200.0){
			psisquared=wf_pp->GetPsiSquared(qmag,rmag,ctheta);
		}
		else
			psisquared=wf_same[id][idprime]->GetPsiSquared(qmag,rmag,ctheta);
	}
	else if(q1q2<0){
		psisquared=wf_opp[id][idprime]->GetPsiSquared(qmag,rmag,ctheta);
	}
	else{
		if(pid==pidprime)
			psisquared=wf_same[id][idprime]->GetPsiSquared(qmag,rmag,ctheta);
		else
			psisquared=1.0;
	}
	
	return psisquared;
}