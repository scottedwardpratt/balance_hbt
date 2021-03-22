#ifndef __BALHBT_H__
#define __BALHBT_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <fstream>
#include "parametermap.h"
#include "resonances.h"
#include "misc.h"
#include "randy.h"
#include "coral.h"
#include "constants.h"
#include <list>

using namespace std;

class CHBTPart;
class Cacceptance;
class CStableInfo;
class CblastWave;
class CHBTCalc;
class CBF;

class CBalHBT{
public:
	int run_number;
	CBalHBT(CparameterMap *parmapin,int run_number);
	CparameterMap *parmap;
	CResList *reslist;
	CblastWave *bw;
	CRandy *randy;
	CBF *bf;
	CHBTCalc *hbtcalc;
	void GetStableInfo(CResList *reslist,double taumax,vector<CStableInfo *> &stablevec,vector<vector<double>> &bfnormweight);
	void GetDecayProducts(CHBTPart *part,vector<CHBTPart *> &products);
	void GetDecayResInfo(CResInfo *resinfo0,int &ndaughters,array<CResInfo *,5> &daughter);
	void Decay(CHBTPart *mother,int &nbodies,array<CHBTPart,5> &daughter);
	void GetPart(vector<CStableInfo *> &stablevec,unsigned int &id);
	void InitHBT(vector<CStableInfo *> &stablevec,string hbtparsfilename);
	void WriteResults(string dirname);
};

class CHBTPart{
public:
	FourVector p;
	FourVector x;
	int motherid;
	CResInfo *resinfo;
	CHBTPart(){};
	~CHBTPart(){};
	CHBTPart(CResInfo *resinfoset){
		resinfo=resinfoset;
	}
	void BjBoost(double dely);
	void Copy(CHBTPart *);
	double GetMass();
	double GetRapidity();
	void Print();
};

class Cacceptance{
public:
	CRandy *randy;
	double ACCEPTANCE;
	Cacceptance(double a){
		ACCEPTANCE=a;
		randy=new CRandy(1234);
	}
	bool Acceptance(CResInfo *resinfo){
		if(randy->ran()<ACCEPTANCE)
			return true;
		else
			return false;
	}
	void Acceptance(CHBTPart *part,bool &acceptQ,bool &acceptP,bool &acceptK,bool &acceptPi,bool &acceptB);
};

class CStableInfo{
public:
	double density;
	CResInfo *resinfo;
	CStableInfo(CResInfo *resinfoset){
		resinfo=resinfoset;
		density=0.0;
	}
	static double denstot;
};
	
class CblastWave{
public:
	CparameterMap *parmap;
	CRandy *randy;
	double Tf,uperpmax,Rperp,etamax,tau;
	double Ybeam,sigmaR;
	double sigma_eta;
	CResList *reslist;
	CblastWave(CparameterMap *parmapin,CRandy *randyset,CResList *reslistset);
	void GenerateParts(vector<CResInfo *> &resinfovec,vector<CHBTPart *> &partvec);
	void GetXP(vector<CHBTPart *> &partvec);
	void GetDecayMomenta(CHBTPart *mother,int &nbodies,vector<CHBTPart *> &daughterpartvec);
	void SetYbeam(double roots);
};

class CHBTCalc{
public:
	CHBTCalc(CparameterMap *parmapin);
	double GetPsiSquared(CHBTPart *part,CHBTPart *partprime);
	CparameterMap *parmap;
	vector<vector<CWaveFunction *>> wf;
};

class CBF{
public:
	CparameterMap *parmap;
	int NYBINS,NPHIBINS;
	double DELPHI,DELY,YMAX;
	vector<double> BFy_pipi;
	vector<double> BFy_piK;
	vector<double> BFy_pip;
	vector<double> BFy_KK;
	vector<double> BFy_Kp;
	vector<double> BFy_pp;
	
	vector<double> BFphi_pipi;
	vector<double> BFphi_piK;
	vector<double> BFphi_pip;
	vector<double> BFphi_KK;
	vector<double> BFphi_Kp;
	vector<double> BFphi_pp;
	
	vector<double> DENOMy_pipi;
	vector<double> DENOMy_piK;
	vector<double> DENOMy_pip;
	vector<double> DENOMy_KK;
	vector<double> DENOMy_Kp;
	vector<double> DENOMy_pp;
	
	vector<double> DENOMphi_pipi;
	vector<double> DENOMphi_piK;
	vector<double> DENOMphi_pip;
	vector<double> DENOMphi_KK;
	vector<double> DENOMphi_Kp;
	vector<double> DENOMphi_pp;
	
	CBF(CparameterMap *parmapin);
	void Zero();
	void Evaluate(vector<CHBTPart *> &partvec,vector<vector<CHBTPart *>> &productvec,
	vector<CHBTPart *> &partprimevec,vector<vector<CHBTPart *>> &productprimevec,double balweight);
	void Increment(CHBTPart *part,CHBTPart *partprime,double weight);
	void PartAntipart(CHBTPart *part);
	void WriteResults(int run_number);
	
	CHBTCalc *hbtcalc;
	CRandy *randy;
};

#endif

