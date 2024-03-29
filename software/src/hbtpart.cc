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
		x[alpha]=part->x[alpha];
	}
	eta=part->eta;
	y=part->y;
	pt=part->pt;
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
	Setp0();
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	y=atanh(p[3]/p[0]);
	eta=asinh(p[3]/pt);
}
void CHBTPart::Setp0(){
	p[0]=sqrt(resinfo->mass*resinfo->mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
}

void CHBTPart::Print(){
	resinfo->Print();
	printf("p=(%g,%g,%g,%g)\n",p[0],p[1],p[2],p[3]);
	printf("x=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3]);
	printf("p^2=%g\n",sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]));
	printf("eta=%g, y=%g, pt=%g\n",eta,y,pt);
}

double CHBTPart::GetDCA(){
	// returns distance of closest approach in cm
	int alpha;
	FourVector y;
	double x2=0.0,pdotx=0.0,p2=0.0,dca=0.0;
	for(alpha=1;alpha<4;alpha++){
		pdotx+=p[alpha]*x[alpha];
		x2+=x[alpha]*x[alpha];
		p2+=p[alpha]*p[alpha];
	}
	for(alpha=0;alpha<4;alpha++){
		y[alpha]=x[alpha]-pdotx*p[alpha]/p2;
		dca+=y[alpha]*y[alpha];
	}
	dca=sqrt(dca)*1.0E-13;
	return dca;
}