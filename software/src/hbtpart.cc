#include "balhbt.h"
using namespace std;

CBalHBT *CHBTPart::balhbt=NULL;

double CHBTPart::GetMass(){
	return sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
}
void CHBTPart::Copy(CHBTPart *part){
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
		p[alpha]=part->p[alpha];
	}
	eta=part->eta;
	y=part->y;
	resinfo=part->resinfo;
}

double CHBTPart::GetRapidity(){
	return atanh(p[3]/p[0]);
}

void CHBTPart::BjBoost(double dely){
	FourVector u;
	u[0]=cosh(dely); u[3]=sinh(dely);
	Misc::Boost(u,p,p);
	Misc::Boost(u,x,x);
}

void CHBTPart::PartAntipart(){
	if(resinfo->q[0]!=0 || resinfo->q[1]!=0 || resinfo->q[2]!=0){
		resinfo=resinfo->reslist->GetResInfoPtr(-resinfo->code);
	}
}

void CHBTPart::SetEtaYPt(){
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	y=atanh(p[3]/p[0]);
	eta=asinh(p[3]/pt);
}

void CHBTPart::Print(){
	resinfo->Print();
	printf("p=(%g,%g,%g,%g)\n",p[0],p[1],p[2],p[3]);
	printf("x=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3]);
	printf("p^2=%g\n",sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]));
	printf("eta=%g, y=%g, pt=%g\n",eta,y,pt);
}