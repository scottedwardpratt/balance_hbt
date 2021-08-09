#ifndef __ACCEPTANCE_BAL_H__
#define __ACCEPTANCE_BAL_H__
using namespace std;
class CHBTPart;
class CBalHBT;

class CAcceptanceBal{
public:
	CAcceptanceBal(CparameterMap *parmap);
	string ACCEPTANCE;
	CBalHBT *balhbt;
	bool acceptance(CHBTPart *part1,CHBTPart *part2,double &efficiency);
	bool acceptance(CHBTPart *part,double &efficiency);
	double etamax,ptmin_pi,ptmin_K,ptmin_p;
};

#endif