#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <cstring>
#include "randy.h"
#include "decay_nbody.h"

using namespace std;

int main(int argc,char *argv[]){
	long long unsigned int im,NM=100000;
	int nbodies=3,iq,NQ=1,ibody;
	printf("Enter nbodies: ");
	scanf("%d",&nbodies);
	printf("Enter NM: ");
	scanf("%lld",&NM);
	CRandy *randy=new CRandy(-time(NULL));
	CDecay_NBody decay(randy);
	vector<FourVector> p(nbodies);
	vector<double> masses(nbodies+1);
	double msum;
	masses[0]=1000.0;
	for (im=0;im<NM;im++){
		msum=0.0;
		for(ibody=1;ibody<=nbodies;ibody++){
			if(randy->ran()<0.5){
				masses[ibody]=0.0;
			}
			else{
				masses[ibody]=randy->ran();
			}
			msum+=masses[ibody];
		}
		double mtarget=masses[0]*randy->ran();
		for(ibody=1;ibody<=nbodies;ibody++){
			masses[ibody]=mtarget*masses[ibody]/msum;
		}
		//sort(masses.rbegin(),masses.rend());
		//masses[0]=0.0;
		//sort(masses.begin(),masses.end());
		//masses[0]=1000.0;
		decay.SetMasses(nbodies,masses);
		for(iq=0;iq<NQ;iq++){
			decay.GenerateMomenta(p);
			/*
			FourVector P={0,0,0,0};
			int alpha;
			for(alpha=0;alpha<4;alpha++)
				P[alpha]=0.0;
			for(ibody=0;ibody<nbodies;ibody++){
				for(alpha=0;alpha<4;alpha++){
					P[alpha]+=p[ibody][alpha];
				}
			}
			printf("P=(%g,%g,%g,%g)\n",P[0],P[1],P[2],P[3]);
			*/
		}
	}
	printf("wmaxmax=%g\n",decay.wmaxmax);
	printf("success rate=%g\n",double(decay.Nsuccess)/double(decay.Ntry));
	return 0;
}
