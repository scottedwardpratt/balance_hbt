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
	CBalHBT *balhbt=new CBalHBT(run_number);
	balhbt->Init();
	printf("Initialization finished\n");
	balhbt->bw->tau=tau;
	balhbt->bw->Rperp=R;
	balhbt->CalcCFs();
	printf("NETWEIGHTsame=%g, NETWEIGHTopp=%g\n",balhbt->bf->NETWEIGHTsame,balhbt->bf->NETWEIGHTopp);
	success=1;
	printf("%d\n",success);
	return 0;
}