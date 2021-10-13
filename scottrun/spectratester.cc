#include <cstdlib>
#include <cstdio>
#include "constants.h"
#include "randy.h"
#include "misc.h"


using namespace std;

int main(int argc,char *argv[]){
	CRandy *randy=new CRandy(-time(NULL));
	double UPERP,T,pt,ptbar=0.0;
	double mass=139.57;
	//double mass=938.27;
	FourVector p,u;
	long long unsigned int imc,NMC=1000000;
	printf("Enter T and UPERP: ");
	scanf("%lf %lf",&T,&UPERP);
	u[2]=u[3]=0.0;
	for(imc=0;imc<NMC;imc++){
		u[1]=UPERP*randy->ran();
		u[2]=UPERP*randy->ran();
		u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]);
		randy->generate_boltzmann(mass,T,p);
		Misc::Boost(u,p);
		pt=sqrt(p[1]*p[1]+p[2]*p[2]);
		ptbar+=pt;
	}
	ptbar=ptbar/double(NMC);
	printf("<pt>=%g\n",ptbar);
	return 0;
}