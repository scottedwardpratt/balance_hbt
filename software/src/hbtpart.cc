#include "balhbt.h"
using namespace std;

double CHBTPart::GetMass(){
	return sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
}
void CHBTPart::Copy(CHBTPart *part){
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
		p[alpha]=part->p[alpha];
	}
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