#include "parametermap.h"
#include "resonances.h"
#include "constants.h"
#include "balhbt.h"

using namespace std;
int main(int argc,char *argv[]){
	int run_number=0;
	CBalHBT *balhbt=new CBalHBT(run_number);
	balhbt->parmap.set("BF_SPECTRA_ONLY",true);
	balhbt->Init();
	balhbt->CalcBFs();
	balhbt->bf->Zero();
	balhbt->CalcBFCFDenoms();
	return 0;
}